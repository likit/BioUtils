'''The script is for splitting sequence in FASTA format to smaller chunks.
Each chunk has an end that overlaps the preceding chunk.

'''

import sys
from pygr import seqdb, sequtil

if len(sys.argv) < 4:
    print >> sys.stderr, \
        'Usage: split_sequence.py fasta_file chunk_size overlap_size'
    raise SystemExit

input_file = sys.argv[1]
chunk_size = int(sys.argv[2])
overlap_size = int(sys.argv[3])

db = seqdb.SequenceFileDB(input_file)

for seq in db:
    window = 0
    if len(db[seq]) <= chunk_size:
        sequtil.write_fasta(sys.stdout, str(db[seq]), id=seq)
    else:
        seq = db[seq]
        _id = 1
        chunk_id = "%s_%d" % (seq.id, _id)
        while window < len(seq):
            chunk = seq[window:window + chunk_size]
            sequtil.write_fasta(sys.stdout, str(chunk), id=chunk_id)
            _id += 1
            window += overlap_size
