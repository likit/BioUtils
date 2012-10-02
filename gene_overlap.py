'''The script reads a gene coordinate from a text file as well as a location
of SNPs or whatever location of interest and report genes that overlap the
location of interest.

'''

from bx.intervals.intersection import Interval, IntervalTree

import sys

def parse_gene_coordinate(infile):
    for line in open(infile):
        cols = line.strip().split(',')
        geneid = cols[0]
        chr, start, end = cols[2:]
        start = int(start)
        end = int(end)
        yield chr, Interval(start, end, value={'geneid':geneid})


def find_overlap(genes, qfile):
    # CHROM  POS     ID      REF     ALT
    for line in open(qfile):
        cols = line.strip().split()
        chrom = cols[0]
        pos = int(cols[1])
        ref = cols[3]
        alt = cols[4]
        insertion = True if len(alt) > len(ref) else False
        if chrom in genes:
            for hit in genes[chrom].find(pos, pos+1):
                print '%s\t%s\t%d\t%d\t%s\t%d\t%s\t%s' % \
                                                    (hit.value['geneid'],
                                                        chrom,
                                                        hit.start,
                                                        hit.end,
                                                        insertion,
                                                        pos,
                                                        ref,
                                                        alt,
                                                    )
                raise SystemExit


def main():
    infile = sys.argv[1]
    qfile = sys.argv[2]

    genes = {}
    for chr, gene in parse_gene_coordinate(infile):
        if chr not in genes:
            genes[chr] = IntervalTree()

        genes[chr].insert_interval(gene)

    find_overlap(genes, qfile)

if __name__=="__main__":
    main()
