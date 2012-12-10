'''The script reads alignments in input_ format and prints out
uninput_sequences sequences in FASTA format.

'''

import sys
from Bio import SeqIO

def parse_input(input_file):
    input_sequences = set()
    for line in open(input_file):
        name = line.strip()
        input_sequences.add(name)
    return input_sequences

def select_seq(input_sequences, fasta_file, min_size):
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id not in input_sequences:
            print >> sys.stderr, 'Writing %s' % record.id
            yield record

def main():
    try:
        input_file = sys.argv[1]
        fasta_file = sys.argv[2]
    except IndexError:
        print >> sys.stderr, \
            'uninput_sequences_seq.py <input_file> <fasta file> [min_size=0]'
        raise SystemExit
    try:
        min_size = sys.argv[3]
    except IndexError:
        min_size = 0

    input_sequences = parse_input(input_file)
    SeqIO.write(select_seq(input_sequences, fasta_file, min_size),
                                                    sys.stdout, 'fasta')


if __name__=='__main__':
    main()
