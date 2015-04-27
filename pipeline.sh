#!/bin/bash
echo `date` ": Starting..."

# MODULES
module load bwa/0.7.3 samtools/0.1.19 python/2.7.3 numpy/1.7.0 scipy/0.12.0b1 pear/0.9.5

## OPTIONS AND DIRECTORIES

READ_1=/net/eichler/vol20/projects/ASD3/backups/fastq/130411_SN711_0331_BD200HACXX/s_5_1.fq.gz
READ_2=/net/eichler/vol20/projects/ASD3/backups/fastq/130411_SN711_0331_BD200HACXX/s_5_2.fq.gz
READ_3=/net/eichler/vol20/projects/ASD3/backups/fastq/130411_SN711_0331_BD200HACXX/s_5_3.fq.gz

BARCODES=/net/eichler/vol20/projects/ASD3/nobackups/BD200HACXX_lane5.txt
BARCODE_LENGTH=8
MOLECULAR_TAG_SIZES="0,5"
MIP_KEY=/net/eichler/vol20/projects/ASD3/nobackups/ASD3_mip_key.txt
OUTPREFIX=SSCP01_SSCP25

FINAL_DEST_DIR=/net/eichler/vol20/projects/ASD3/nobackups/nkrumm

MIPGEN=/net/eichler/vol20/projects/ASD3/nobackups/nkrumm/MIPGEN/
BWA_INDEX_SOURCE=/net/eichler/vol3/home/bcoe/hg19/hg19_masked_arhgap11b_numeric.fa
THREADS=6

# TEMP DIR SETUP
TMP=$(mktemp -d)

# COPY SOURCE to $TMP
rsync -arv --bwlimit=15000 $BWA_INDEX_SOURCE* $TMP


echo `date` ": Copying read FASTQ data to local temp directory $TMP)"

rsync $READ_1 $TMP
rsync $READ_2 $TMP
rsync $READ_3 $TMP

READ_1_TMP=$TMP/`basename $READ_1`
READ_2_TMP=$TMP/`basename $READ_2`
READ_3_TMP=$TMP/`basename $READ_3`

# SCRIPT
echo `date` ": Running mipgen_fq_cutter_pe.py..."
python $MIPGEN/tools/mipgen_fq_cutter_pe.py \
    $READ_1_TMP $READ_3_TMP -i $READ_2_TMP \
    -j $BARCODE_LENGTH \
    -tb $BARCODES \
    -o $TMP/$OUTPREFIX.barcoded

echo `date` ": Running pear..."
pear -f $TMP/$OUTPREFIX.barcoded.r1.indexed.fq \
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
echo `date` ": Saving results to $FINAL_DEST_DIR"
rsync -arv $TMP/$OUTPREFIX.barcoded.indexed.sort.bam $FINAL_DEST_DIR

echo `date` ": Done."
