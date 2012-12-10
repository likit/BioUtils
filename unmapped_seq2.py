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
        print >> sys.stderr, 'unmapped_seq2.py <text file> <fasta file>'

    db = seqdb.SequenceFileDB(fasta_file)

    input_sequences = set()
    print >> sys.stderr, 'Reading sequences...'
    for line in open(input_file):
        input_sequences.add(line.strip())

    print >> sys.stderr, 'Writing unmapped sequences...'
    for seq in db:
        if seq not in input_sequences:
            sequtil.write_fasta(sys.stdout, str(db[seq]), id=seq)

if __name__=='__main__':
    main()
