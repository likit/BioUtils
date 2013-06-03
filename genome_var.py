'''The script generates edited copy of a genome by replacing
reference bases with SNPs from VCF files.

'''

import sys
import csv

from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord

def parse_vcf(varfile):
    reader = csv.reader(open(varfile), "excel-tab")
    for line in reader:
        if line[0][0] == "#":
            continue

        pos = int(line[1]) - 1
        var = line[4].split(',')

        yield pos, var

for seq_record in SeqIO.parse(sys.argv[1], 'fasta'):
    print >> sys.stderr, "Seq ID = %s, Length = %d" % \
                                        (seq_record.id, len(seq_record))
    seq = MutableSeq(str(seq_record.seq))

    n = 0
    for pos, var in parse_vcf(sys.argv[2]):
        # if (len(var) > 2) or (len(var[0]) > 1):
            # continue
        if (len(var) > 1) or (len(var[0]) > 1):
            continue
        else:
            seq[pos] = var[0]
            n += 1

    SeqIO.write(SeqRecord(Seq(str(seq)), id=seq_record.id),
                                            sys.stdout, 'fasta')

    print >> sys.stderr, "Total variants = %d" % n
    print >> sys.stderr, "Percent variation = %.3f" % (float(n)*100/len(seq))

