'''Select sequences from FASTA file according to a given list
in a text file or standard input.

Output is written to standard output.

'''

import sys
from Bio import SeqIO


def getids(input_file):
    thelist = []
    for seqid in input_file:
        thelist.append(seqid.strip())
    return thelist

def select(thelist, fasta_file):
    for rec in SeqIO.parse(fasta_file, 'fasta'):
        if rec.id in thelist:
            SeqIO.write(rec, sys.stdout, 'fasta')

if __name__=='__main__':
    try:
        input_file = open(sys.argv[1])
        fasta_file = sys.argv[2]
    except IndexError, IOError:
        print >> sys.stderr, 'Usage: python get_fasta_seq.py' + \
                                '<list file> <fasta file>'
        raise SystemExit
    else:
        thelist = getids(input_file)
        select(thelist, fasta_file)
