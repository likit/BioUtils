'''The script reads sequences in fasta format and selects sequences that are
longer than a user specified length.

'''

import sys
from Bio import SeqIO

def main():
    try:
        fasta_file = sys.argv[1]
        size = int(sys.argv[2])
    except IndexError:
        print >> sys.stderr, 'size_select.py <fasta file> <size>'

    for record in SeqIO.parse(fasta_file, 'fasta'):
        if len(record) >= size:
            print >> sys.stderr, 'Writing %s ..' % (record.id)
            SeqIO.write(record, sys.stdout, 'fasta')


if __name__=='__main__':
    main()
