'''The script reads alignments in PSL format and prints out
unmapped sequences in FASTA format.

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
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id not in mapped_seqs:
            print >> sys.stderr, record.id
            yield record

if __name__=='__main__':
    try:
        pslfile = sys.argv[1]
        fastafile = sys.argv[2]
    except IndexError:
        print >> sys.stderr, 'unmapped_seq.py <psl file> <fasta file>'
    mapped_seqs = parse_mapped(pslfile)
    print >> sys.stderr, 'Writing unmapped sequences ...'
    SeqIO.write(select_seq(mapped_seqs, fastafile), sys.stdout, 'fasta')
