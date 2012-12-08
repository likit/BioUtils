'''The script reads results from BLAST alignments in blast9 format and
store transcript ID from Ensembl with a matching ID from given gene models
in BED format.

'''

import sys
import cPickle

if len(sys.argv) < 3:
    print >> sys.stderr, "Usage: gene_model_match.py blast9_file output_file"

trans_db = {}  # transcripts database

print >> sys.stderr, "Reading %s.." % sys.argv[1]

for line in open(sys.argv[1]):
    if line.startswith("#"): continue
    row = line.strip().split()
    tid = row[0]
    ensid = row[1]
    if tid in trans_db:
        continue
    else:
        trans_db[tid] = ensid

print >> sys.stderr, "Writing transcripts database to %s" % (sys.argv[2])
cPickle.dump(trans_db, open(sys.argv[2], "w"))
print >> sys.stderr, "Done."
print >> sys.stderr, "Total transcripts = %d" % (len(trans_db))
