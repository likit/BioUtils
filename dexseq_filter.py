'''The script reads results from DEXSeq pacakge and filters out DEU located
at either 3' or 5' end of any gene model.

How to obtain coord.txt file from R.

> library("DEXSeq")
> coord <- fData(expData)[,c(1,2,9:12)]
> write.table(coord, file="coord.txt", sep="\t",
                    row.names=F, col.names=F, quote=F)

where expData is ExonCountSet object.

'''

import csv
import sys
import itertools

from collections import namedtuple
from bx.intervals.intersection import Interval, IntervalTree

ExonQuery = namedtuple('ExonQuery', ['chrom', 'start', 'end'])

def parse_models(bedfile):
    '''Converts gene models in BED format to a list of exon intervals.'''

    reader = csv.reader(open(bedfile), dialect='excel-tab')
    for row in reader:
        exons = []
        chrom = row[0]
        chrom_start = int(row[1])
        geneid = row[3]

        exon_sizes = [int(s) for s in row[10].split(',')]
        exon_starts = [chrom_start + int(s) for s in row[11].split(',')]

        for i in range(len(exon_starts)):
            if i == 0:
                terminal = True
            elif i == len(exon_starts) - 1:
                terminal = True
            else:
                terminal = False

            exon_start = exon_starts[i]
            exon_end = exon_start + exon_sizes[i]

            exon = Interval(exon_start, exon_end, value={'geneid':geneid,
                                                'terminal':terminal,
                                                'chrom':chrom})
            exons.append(exon)

        yield exons

def is_terminal(genes, exon):
    for hit in genes[exon.chrom].find(exon.start, exon.end):
        if hit.value['terminal']: return True

    return False

def add_intervals(genes, exons):
    '''Add exons to interval tree.'''
    for exon in exons:
        if exon.value['chrom'] not in genes:
            genes[exon.value['chrom']] = IntervalTree()

        genes[exon.value['chrom']].insert_interval(exon)

def parse_dexseq(result_file, coord_file):
    '''Reads results from DEXSeq package and return ExonQuery.'''

    results = csv.reader(open(result_file), dialect='excel-tab')
    coords = csv.reader(open(coord_file), dialect='excel-tab')
    for res_row, coord_row in itertools.izip(results, coords):
        chrom = coord_row[2]
        start = int(coord_row[3])
        end = int(coord_row[4])
        yield ExonQuery(chrom, start, end), res_row

def main():
    if len(sys.argv) < 4:
        print >> sys.stderr, 'Usage: dexseq_filter.py bed_file ' + \
                                            'result_file coord_file'
        raise SystemExit

    bedfile = sys.argv[1]
    result_file = sys.argv[2]
    coord_file = sys.argv[3]
    genes = {}
    for n, exons in enumerate(parse_models(bedfile), start=1):
        add_intervals(genes, exons)

    # print is_terminal(genes, ExonQuery('chr10', 3160740, 3160899))

    for exon_query, result_row in parse_dexseq(result_file, coord_file):
        if is_terminal(genes, exon_query):
            result_row.append('YES')
        else:
            result_row.append('NO')

        print '\t'.join(result_row)

if __name__=='__main__':
    main()
