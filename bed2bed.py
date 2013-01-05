'''Convert gene models in BED format to BED format to be used as an
input of multiBamCov to count reads map to each transcripts when reads are
mapped directly to transcript sequences not genome sequence.
Output is written to standard output.

'''

import sys
import csv

input_file = sys.argv[1]
reader = csv.reader(open(input_file), dialect='excel-tab')
for row in reader:
    gene_id = row[3]
    size = sum([int(s) for s in row[10].split(',')])
    strand = row[5]
    writer = csv.writer(sys.stdout, dialect='excel-tab')
    writer.writerow([gene_id, 0, size, gene_id, 1000,
                        strand, 0, size, '0,0,0', 1, size, 0])
