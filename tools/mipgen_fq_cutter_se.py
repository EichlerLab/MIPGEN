# Written by Evan Boyle
# boylee [at] uw.edu
import sys
import re
from subprocess import Popen, PIPE
from optparse import OptionParser
from distutils import spawn

# for MacOS compatability
ZCAT_BINARY=spawn.find_executable("zcat") or spawn.find_executable("gzcat")

def fastq_reader(filename):
    """Yields a read from a fastq file, optionally gzipped"""
    if options.gzip_input or re.search(".gz$", filename):
        f = Popen(ZCAT_BINARY + ' ' + filename, shell=True, stdout=PIPE).stdout
    else:
        f = open(filename, 'r')
    for line in f:
        qname = line.rstrip("\n")
        seq = f.next().rstrip("\n")
        sep = f.next().rstrip("\n")
        qual = f.next().rstrip("\n")
        yield [qname, seq, sep, qual]

parser = OptionParser("%prog read.fq [options]")
parser.add_option("-o", "--output_prefix", type="str",
                  help="directs output to given path")
parser.add_option("-i", "--index_read", dest="index_file", type="str",
                  help="index read file")
parser.add_option("-j", "--index_length", type="int",
                  help="truncate index read to given length")
parser.add_option("-b", "--barcode_file", type="str",
                  help="select barcodes in file of label<tab>sequence")
parser.add_option("-t", "--tolerant", action="store_true", default=False,
                  help="allow 1bp edits from provided barcodes in file")
parser.add_option("-m", "--molecular_tag", type="str",
                  help="molecular tag, template: \"-m 3,2\" for first 3 bases, last 2 bases")
parser.add_option("-l", "--truncated_read_length", type="int",
                  help="truncate reads to given length")
parser.add_option("-L", "--skip_length", type="int",
                  help="skip the number of bases provided (post truncation)")
parser.add_option("-d", "--discard", action="store_true",  default=False,
                  help="discard read pairs for which one of the reads is all '#' quality")
parser.add_option("-e", "--partial_discard", action="store_true", default=False,
                  help="enables truncation of reads to eliminate tailing bases of '#' quality")
parser.add_option("-z", "--gzip_input", action="store_true", default=False,
                  help="input fastq files are gzipped; done by default if file ends in \".gz\"")

(options, args) = parser.parse_args()

if(options.partial_discard and not options.discard):
    print "-e option requires -d option"
if(options.index_file is None and options.index_length is not None):
    print "-j option requires -i option"

used_barcodes = set()
bases = "ATCGN"

if options.barcode_file:
    with open(options.barcode_file, 'r') as b_in:
        for barcode_line in b_in:
            (barcode_label, barcode_seq) = barcode_line.rstrip().split()
            used_barcodes.add(barcode_seq)
            if options.tolerant:
                for i in range(len(barcode_seq)):
                    native_seq = list(barcode_seq)
                    for j in range(5):
                        mutated_seq = native_seq
                        mutated_seq[i] = bases[j]
                        used_barcodes.add("".join(mutated_seq))

if options.molecular_tag is not None:
    molecular_tag_specs = options.molecular_tag.split(",")
    molecular_tag_specs = [int(entry) for entry in molecular_tag_specs]

first_in = fastq_reader(sys.argv[1])

if(options.output_prefix is not None):
    outfq = options.output_prefix + ".indexed.fq"
else:
    outfq = sys.argv[1] + ".indexed.fq"

if options.index_file:
    i_in = fastq_reader(options.index_file)

with open(outfq, 'w') as out:
    while(1):
        try:
            block = first_in.next()
        except StopIteration:
            # end of file
            break

        block[0] = re.sub("/\d$", "", block[0])

        if options.index_file:
            index_block = i_in.next()
            barcode = index_block[1]
        else:
            barcode_in_header = re.search("#([ATGCN]+)(-[ATGCN]+)?$", block[0])
            if barcode_in_header:
                barcode = barcode_in_header.group(1)
            else:
                barcode = "N"

        if options.index_file and options.index_length:
            barcode = barcode[:options.index_length]

        if options.barcode_file and barcode not in used_barcodes:
            continue

        if options.index_file is None and barcode_in_header:
            pass
        else:
            block[0] += "#" + barcode
            block[0] = block[0].replace(" ", "_")
            if options.molecular_tag:
                tag_in_header = re.search("#[ATGCN]+-[ATGCN]+$", block[0])

        if options.molecular_tag and not tag_in_header:
            tag = block[1][:molecular_tag_specs[0]] \
                + block[1][len(block[1]) - molecular_tag_specs[1] - 1: -1]
            block[0] += "-" + tag
            block[1] = block[1][molecular_tag_specs[0]:len(block[1]) - 1 - molecular_tag_specs[1]]
            block[3] = block[3][molecular_tag_specs[0]:len(block[3]) - 1 - molecular_tag_specs[1]]

        if options.truncated_read_length:
            block[1] = block[1][:options.truncated_read_length]
            block[3] = block[3][:options.truncated_read_length]

        if options.partial_discard:
            match = re.search("#+$", block[3])
            if match:
                block[1] = block[1][0:match.start()]
                block[3] = block[3][0:match.start()]

        if options.skip_length is not None:
            block[1] = block[1][options.skip_length:]
            block[3] = block[3][options.skip_length:]

        out.write("%s\n%s\n+\n%s\n" % (block[0], block[1], block[3]))

first_in.close()

if options.index_file is not None:
    i_in.close()

if options.output_prefix is not None:
    samfile = options.output_prefix + ".indexed.sam"
else:
    suffix_match = re.search("fq$", sys.argv[1])
    if suffix_match is None:
        suffix_match = re.search("fastq$", sys.argv[1])
    if suffix_match is not None:
        samfile = sys.argv[1][:suffix_match.start()] + "indexed.sam"
    else:
        samfile = sys.argv[1] + ".indexed.sam"

sys.stderr.write("#fq cutting finished\n")
sys.stderr.write("#suggested output file for alignment: " + samfile + "\n")
