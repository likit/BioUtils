'''Extract only alignments of sequences in a list file.'''

import sys

pslfile = sys.argv[1]
listfile = sys.argv[2]
try:
    seqtype = sys.argv[3]
except IndexError:
    seqtype = 9
else:
    if seqtype == 'query':
        seqtype = 9
    elif seqtype == 'target':
        seqtype = 13
    else:
        print >> sys.stderr, 'Unregconized sequence type.'
        raise SystemExit

sequences = set([seq.strip() for seq in open(listfile)])

for align in open(pslfile):
    query = align.split()[seqtype]
    if query in sequences:
        print align,
