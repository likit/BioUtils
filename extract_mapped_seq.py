'''Reads alignments in PSL format and FASTA format and extract a region or
a sequence that do not mapped to a reference sequence.
A script requires a -singleHit parameter for pslRep.

'''

import sys
from collections import namedtuple
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

SeqAlign = namedtuple('SeqAlign', ['name', 'chrom', 'qStart', 'qEnd'])

def parse_psl(psl_file):
    print >> sys.stderr, 'parsing %s...' % psl_file
    alignments = {}
    for n, line in enumerate(open(psl_file)):
        cols = line.strip().split()
        name = cols[9]
        qStart = int(cols[11])
        qEnd = int(cols[12])
        chrom = cols[13]
        seqalign = SeqAlign(name, chrom, qStart, qEnd)
        alignments[name] = seqalign

        if n % 1000 == 0:
            print >> sys.stderr, '...', n

    return alignments

def extract(alignments, fasta_file, mincov=30, minlen=100):
    print >> sys.stderr, 'parsing %s...' % fasta_file
    for n, rec in enumerate(SeqIO.parse(fasta_file, 'fasta')):
        seq = rec.seq.tostring()
        try:
            algn = alignments[rec.name]
        except KeyError: # an entire sequence is unmapped
            SeqIO.write(rec, sys.stdout, 'fasta')
        else:
            mapped_seq = seq[algn.qStart:algn.qEnd]

            if float(len(mapped_seq))/len(rec) * 100 < mincov:
                continue

            head, sep, tail = seq.partition(mapped_seq)
            if head and len(head) >= minlen:
                seqrec = SeqRecord(id=rec.name+'_1',
                                    description='',
                                    seq=Seq(head))
                SeqIO.write(seqrec, sys.stdout, 'fasta')
            if tail and len(tail) >= minlen:
                seqrec = SeqRecord(id=rec.name+'_2',
                                    description='',
                                    seq=Seq(tail))
                SeqIO.write(seqrec, sys.stdout, 'fasta')

        if n % 1000 == 0:
            print >> sys.stderr, '...', n

if __name__=='__main__':
    if len(sys.argv) < 5:
        print >> sys.stderr, 'Usage: python extract_unmapped_seq.py' + \
                ' <psl file> <fasta file> [min coverage]' + \
                '[min unmapped length]'

        raise SystemExit
    else:
        psl_file = sys.argv[1]
        fasta_file = sys.argv[2]
        mincov = int(sys.argv[3])
        minlen = int(sys.argv[4])

        alignments = parse_psl(psl_file)
        extract(alignments, fasta_file, mincov, minlen)
