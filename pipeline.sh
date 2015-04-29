#!/bin/bash
set -e

echo `date` ": Starting..."

## OPTIONS AND DIRECTORIES

# FASTQ_FOLDER=$1
# FLOWCELL=$2
# LANE=$3
# BARCODES=$4
OUTPREFIX=${FLOWCELL}_${LANE}

READ_1=${FASTQ_FOLDER}/s_${LANE}_1.fq.gz
READ_2=${FASTQ_FOLDER}/s_${LANE}_2.fq.gz
READ_3=${FASTQ_FOLDER}/s_${LANE}_3.fq.gz


BARCODE_LENGTH=8
MOLECULAR_TAG_SIZES="0,5"
MIP_KEY=/net/eichler/vol20/projects/ASD3/nobackups/ASD3_mip_key.txt
LOCK_FILE=/net/eichler/vol20/projects/ASD3/nobackups/_lockfile
FINAL_DEST_DIR=/net/eichler/vol20/projects/ASD3/nobackups/analysis_complete_nkrumm

MIPGEN=/net/eichler/vol20/projects/ASD3/nobackups/nkrumm/MIPGEN/
BWA_INDEX_SOURCE=/net/eichler/vol20/projects/ASD3/nobackups/nkrumm/hg19_index/hg19_masked_arhgap11b_numeric.fa
THREADS=6

# TEMP DIR SETUP
TMP=$(mktemp -d)

LOCK_COUNT=$(cat $LOCK_FILE)

while [[ LOCK_COUNT > 4 ]]; do
    echo "Waiting for rsync lock;"
    sleep $(( ( RANDOM % 100 )  + 60 ))  # wait a randomish amount of time
    LOCK_COUNT=$(cat $LOCK_FILE)
done

LOCK_COUNT=$((LOCK_COUNT+1))
echo $LOCK_COUNT > $LOCK_FILE

# COPY SOURCE to $TMP
rsync -arv $BWA_INDEX_SOURCE* $TMP


echo `date` ": Copying read FASTQ data to local temp directory $TMP)"

rsync $READ_1 $TMP
rsync $READ_2 $TMP
rsync $READ_3 $TMP

READ_1_TMP=$TMP/`basename $READ_1`
READ_2_TMP=$TMP/`basename $READ_2`
READ_3_TMP=$TMP/`basename $READ_3`

LOCK_COUNT=$(cat $LOCK_FILE)
LOCK_COUNT=$((LOCK_COUNT-1))
echo $LOCK_COUNT > $LOCK_FILE


# MODULES
module load bwa/0.7.3 samtools/0.1.19 python/2.7.3 numpy/1.7.0 scipy/0.12.0b1 pear/0.9.5

# SCRIPT

echo `date` ": Running mipgen_fq_cutter_pe.py..."
python $MIPGEN/tools/mipgen_fq_cutter_pe.py \
    $READ_1_TMP $READ_3_TMP -i $READ_2_TMP \
    -j $BARCODE_LENGTH \
    -tb $BARCODES \
    -o $TMP/$OUTPREFIX.barcoded

echo `date` ": Running pear..."
pear --threads $THREADS --memory 4G \
    -f $TMP/$OUTPREFIX.barcoded.r1.indexed.fq \
    -r $TMP/$OUTPREFIX.barcoded.r2.indexed.fq \
    -o $TMP/$OUTPREFIX.barcoded

echo `date` ": Running mipgen_fq_cutter_se.py..."
python $MIPGEN/tools/mipgen_fq_cutter_se.py \
    $TMP/$OUTPREFIX.barcoded.assembled.fastq \
    -tb $BARCODES \
    -m $MOLECULAR_TAG_SIZES \
    -o $TMP/$OUTPREFIX.barcoded

echo `date` ": Running bwa and sorting..."
bwa mem -t $THREADS $TMP/`basename $BWA_INDEX_SOURCE` \
    $TMP/$OUTPREFIX.barcoded.indexed.fq \
    | samtools view -bS - \
    | samtools sort - $TMP/$OUTPREFIX.barcoded.indexed.sort

echo `date` ": Running mipgen_smmip_collapser.py..."
samtools view -h $TMP/$OUTPREFIX.barcoded.indexed.sort.bam \
    | python $MIPGEN/tools/mipgen_smmip_collapser.py 5 \
        $TMP/$OUTPREFIX.barcoded.indexed.sort.collapse \
        -m $MIP_KEY \
        -s -f 1 -r -b \
        $BARCODES

# COPY OUT FINAL
echo `date` ": Saving results to $FINAL_DEST_DIR/$OUTPREFIX"
mkdir -p $FINAL_DEST_DIR/$OUTPREFIX
rsync -arv $TMP/$OUTPREFIX.barcoded.indexed.sort.* $FINAL_DEST_DIR/$OUTPREFIX/

# cleaning up
rm -rf $TMP
echo `date` ": Done."
