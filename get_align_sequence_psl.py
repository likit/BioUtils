'''Get a part of a target sequence that aligned to a given query sequence.'''

import sys
import csv
from pygr import seqdb, sequtil

psl_file = sys.argv[1]
genome_file = sys.argv[2]
genome = seqdb.SequenceFileDB(genome_file)

reader = csv.reader(open(psl_file))
for cols in reader:
    target = cols[13]
    start = int(cols[15])
    end = int(cols[16])
    seq = genome[target][start:end]
    seqid = target + '_' + cols[9]
    sequtil.write(seq, sys.stdout, id=seqid)
