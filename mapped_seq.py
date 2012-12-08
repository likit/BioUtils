'''The script reads alignment results from BLAT in PSL format and print out
sequences of those sequences in FASTA format.

'''

import sys
from Bio import SeqIO

def parse_mapped(pslfile):
    mapped = set()
    for line in open(pslfile):
        name = line.strip().split()[9]
        mapped.add(name)
    return mapped

def select_seq(mapped_seqs, fasta_file):
    for n, record in enumerate(SeqIO.parse(fasta_file, 'fasta'), start=1):
        if record.id in mapped_seqs:
            yield record
        if n % 1000 == 0:
            print >> sys.stderr, '...', n

if __name__=='__main__':
    try:
        pslfile = sys.argv[1]
        fastafile = sys.argv[2]
    except IndexError:
        print >> sys.stderr, 'unmapped_seq.py <psl file> <fasta file>'
    mapped_seqs = parse_mapped(pslfile)
    print >> sys.stderr, 'Writing mapped sequences ...'
    SeqIO.write(select_seq(mapped_seqs, fastafile), sys.stdout, 'fasta')
