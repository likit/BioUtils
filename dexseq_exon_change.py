'''The script reads results from DEXSeq and append a column with the number
of DEU for each genes with padjust < a given FDR.

'''

import sys

if len(sys.argv) == 1:
    print >> sys.stderr, 'Usage: dexseq_exon_change.py <result_file> [FDR]'
    print >> sys.stderr, 'Default FDR = 0.1'
    raise SystemExit

result_file = sys.argv[1]
try:
    FDR = float(sys.argv[2])
except IndexError:
    FDR = 0.1  # default value

genes = {}
for line in open(result_file):
    row = line.strip().split()
    name = row[0]

    if row[4] == "NA":
        genes[name] = genes.get(name, 0)
        continue

    if float(row[4]) < FDR:
        genes[name] = genes.get(name, 0) + 1
    else:
        genes[name] = genes.get(name, 0)

for line in open(result_file):
    row = line.strip().split()
    name = row[0]
    row.append(str(genes[name]))
    print '\t'.join(row)
