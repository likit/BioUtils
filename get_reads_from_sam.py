#! /usr/local/bin/python

'''Reads BAM file and filter out unmapped reads and write them in FASTA
format. Paired-end reads without proper mate-pair are filtered out.
Reads are written to standard output.

Author: Likit Preeyanon
Email: preeyano@msu.edu

'''

import pysam
import sys

def write_fasta(read):
    print '>%s\n%s' % (read.qname, read.seq)

def writeReads(infile):
    samfile = pysam.Samfile(infile, 'rb')
    for read in samfile.fetch():
        if read.is_unmapped:
            continue
        if read.is_paired:
            if read.is_proper_pair:
                if read.is_read1:
                    read.qname = read.qname + '/1'
                elif read.is_read2:
                    read.qname = read.qname + '/2'
                else:
                    raise ValueError, 'Unrecognized read'
                write_fasta(read)
        else:
            write_fasta(read)


if __name__=='__main__':

    try:
        writeReads(sys.argv[1])
    except IOError, IndexError:
        print >> sys.stderr, 'Usage: python get_reads_from_sam.py <bam file>'
        raise SystemExit
