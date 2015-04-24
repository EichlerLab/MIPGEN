# Written by Evan Boyle
# boylee [at] uw.edu

import sys
import re
from subprocess import Popen, PIPE
from optparse import OptionParser


def fastq_reader(filename):
    """Yields a read from a fastq file, optionally gzipped"""
    if options.gzip_input or re.search(".gz$", filename):
        f = Popen('gzcat ' + filename, shell=True, stdout=PIPE).stdout
    else:
        f = open(filename, 'r')
    for line in f:
        qname = line.rstrip("\n")
        seq = f.next().rstrip("\n")
        sep = f.next().rstrip("\n")
        qual = f.next().rstrip("\n")
        yield [qname, seq, sep, qual]

parser = OptionParser("%prog read_1.fq read_2.fq [options]")
parser.add_option("-o", "--output_prefix", dest="output_prefix", type="str",
                  help="directs output to given path")
parser.add_option("-i", "--index_read", dest="index_file", type="str",
                  help="index read file to insert barcode into header")
parser.add_option("-j", "--index_length", dest="index_length", type="int",
                  help="truncate index read to given length")
parser.add_option("-b", "--barcode_file", type="str",
                  help="select barcodes in file of label<tab>sequence")
parser.add_option("-t", "--tolerant", action="store_true", default=False,
                  help="allow 1bp edits from provided barcodes in file")
parser.add_option("-m", "--molecular_tag", dest="molecular_tag", type="str",
                  help="molecular tag, template: \"-m 10,0\" for read 1, 10bp")
parser.add_option("-l", "--truncate_read_length", dest="truncated_read_length", type="int",
                  help="truncate reads to given length")
parser.add_option("-L", "--skip_length", dest="skip_length", type="int",
                  help="skip the number of bases provided (post truncation)")
parser.add_option("-d", "--discard", action="store_true", dest="discard", default=False,
                  help="discard read pairs for which one of the reads is all '#' quality")
parser.add_option("-e", "--partial_discard", action="store_true", dest="partial_discard",
                  default=False, help="enables truncation of reads to eliminate tailing bases of '#' quality")  # noqa
parser.add_option("-z", "--gzip_input", action="store_true", dest="gzip_input", default=False,
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

if options.molecular_tag:
    molecular_tag_specs = options.molecular_tag.split(",")
    molecular_tag_specs = [int(entry) for entry in molecular_tag_specs]

# open files
first_in = fastq_reader(sys.argv[1])
second_in = fastq_reader(sys.argv[2])
if options.index_file:
    i_in = fastq_reader(options.index_file)

if options.output_prefix:
    first_outfq = options.output_prefix + ".r1.indexed.fq"
    second_outfq = options.output_prefix + ".r2.indexed.fq"
else:
    first_outfq = sys.argv[1] + ".r1.indexed.fq"
    second_outfq = sys.argv[2] + ".r2.indexed.fq"

# MAIN LOOP

with open(first_outfq, 'w') as first_out:
    with open(second_outfq, 'w') as second_out:
        while(1):
            try:
                first_block = first_in.next()
                second_block = second_in.next()
            except StopIteration:
                # end of file
                break
            if (options.discard and (re.match("^#+$", first_block[3].rstrip())
                                     or re.match("^#+$", second_block[3].rstrip()))):
                continue
            first_block[0] = re.sub("/\d$", "", first_block[0], 1)
            second_block[0] = re.sub("/\d$", "", second_block[0], 1)

            if len(first_block[0]) >= len(second_block[0]):
                # header must be same for each read; take the longer
                second_block[0] = first_block[0]
            else:
                first_block[0] = second_block[0]

            if options.index_file:
                index_block = i_in.next()
                barcode = index_block[1]
            else:
                barcode_in_first_header = re.search("#([ATGCN]+)(-[ATGCN]*)?$", first_block[0])
                barcode_in_second_header = re.search("#([ATGCN]+)(-[ATGCN]*)?$", second_block[0])
                if (barcode_in_first_header and barcode_in_second_header
                   and barcode_in_first_header.group(1) == barcode_in_second_header.group(1)):
                    barcode = barcode_in_first_header.group(1)
                else:
                    barcode = "N"

            if (options.index_file and options.index_length):
                barcode = barcode[:options.index_length]

            if (options.barcode_file and barcode not in used_barcodes):
                continue

            if (options.index_file is None and barcode_in_first_header
                    and barcode_in_second_header):
                pass
            else:
                first_block[0] += "#" + barcode
                first_block[0] = first_block[0].replace(" ", "_")
                second_block[0] += "#" + barcode
                second_block[0] = second_block[0].replace(" ", "_")
                if options.molecular_tag:
                    tag_in_first_header = re.search("#[ATGCN]+-[ATGCN]+$", first_block[0])
                    tag_in_second_header = re.search("#[ATGCN]+-[ATGCN]+$", second_block[0])

            if options.molecular_tag:
                if not tag_in_first_header and not tag_in_second_header:
                    tag = first_block[1][:molecular_tag_specs[0]] \
                        + second_block[1][:molecular_tag_specs[1]]
                    first_block[0] += "-" + tag
                    first_block[1] = first_block[1][molecular_tag_specs[0]:]
                    first_block[3] = first_block[3][molecular_tag_specs[0]:]
                    second_block[0] += "-" + tag
                    second_block[1] = second_block[1][molecular_tag_specs[1]:]
                    second_block[3] = second_block[3][molecular_tag_specs[1]:]
                elif not (tag_in_first_header and tag_in_second_header):
                    sys.stderr.write("molecular tag formatting anomaly; canceling analysis\n")
                    break
                else:  # tag already present
                    pass  # do nothing

            if options.truncated_read_length:
                first_block[1] = first_block[1][:options.truncated_read_length]
                first_block[3] = first_block[3][:options.truncated_read_length]
                second_block[1] = second_block[1][:options.truncated_read_length]
                second_block[3] = second_block[3][:options.truncated_read_length]

            if options.partial_discard:
                match = re.search("#+$", first_block[3])
                if match:
                    first_block[1] = first_block[1][0:match.start()]
                    first_block[3] = first_block[3][0:match.start()]
                match = re.search("#+$", second_block[3])
                if match:
                    second_block[1] = second_block[1][0:match.start()]
                    second_block[3] = second_block[3][0:match.start()]

            if options.skip_length:
                first_block[1] = first_block[1][options.skip_length:]
                first_block[3] = first_block[3][options.skip_length:]
                second_block[1] = second_block[1][options.skip_length:]
                second_block[3] = second_block[3][options.skip_length:]

            # write output
            first_out.write("%s\n%s\n+\n%s\n" %
                            (first_block[0], first_block[1], first_block[3]))
            second_out.write("%s\n%s\n+\n%s\n" %
                             (second_block[0], second_block[1], second_block[3]))

first_in.close()
second_in.close()

if options.index_file:
    i_in.close()

if options.output_prefix:
    samfile = options.output_prefix + ".indexed.sam"
else:
    f_shared_index = 0

    if(sys.argv[1][-1] == sys.argv[2][-1]):
        r_shared_index = -1
    else:
        r_shared_index = None

    while(sys.argv[1][f_shared_index] == sys.argv[2][f_shared_index]):
        f_shared_index += 1

    while(sys.argv[1][r_shared_index - 1] == sys.argv[2][r_shared_index - 1]):
        r_shared_index -= 1
    root_component = sys.argv[1][:f_shared_index]
    derived_component = sys.argv[1][f_shared_index:r_shared_index] \
        + sys.argv[2][f_shared_index:r_shared_index]
    samfile = root_component + derived_component + ".indexed.sam"

sys.stderr.write("#fq cutting finished\n")
sys.stderr.write("#suggested output file for alignment: " + samfile + "\n")
