'''The script reads a list of sequences and prints out sequences not in a
list in FASTA format.

'''

import sys
from pygr import seqdb, sequtil

def main():
    try:
        input_file = sys.argv[1]
        fasta_file = sys.argv[2]
    except IndexError:
        print >> sys.stderr, \
            'unmapped_seq2.py <text file> <fasta file> [min length=50]'
    try:
        min_length = int(sys.argv[3])
    except IndexError:
        min_length = 50

    db = seqdb.SequenceFileDB(fasta_file)

    input_sequences = set()
    print >> sys.stderr, 'Reading sequences...',
    for line in open(input_file):
        input_sequences.add(line.strip())
    print >> sys.stderr, 'total %d' % len(input_sequences)

    print >> sys.stderr, 'Writing unmapped sequences...'
    for seq in db:
        sequence = db[seq]
        if (seq not in input_sequences and len(sequence) > min_length):
            sequtil.write_fasta(sys.stdout, str(sequence), id=seq)
            print >> sys.stderr, '%s has been written' % seq


if __name__=='__main__':
    main()
