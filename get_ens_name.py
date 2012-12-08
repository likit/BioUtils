'''The script reads a list of transcripts from standard input and reports
Ensembl ID from pickle transcripts database file.

'''

import sys
import cPickle
import csv

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
            exon_start = exon_starts[i]
            exon_end = exon_start + exon_sizes[i]

            exon = Interval(exon_start, exon_end, value={
                                                        'geneid':geneid,
                                                        'chrom':chrom},
                                                        )
            exons.append(exon)
        yield exons

def add_intervals(genes, exons):
    '''Add exons to interval tree.'''
    for exon in exons:
        if exon.value['chrom'] not in genes:
            genes[exon.value['chrom']] = IntervalTree()

        genes[exon.value['chrom']].insert_interval(exon)

def find_id(genes, query, trans_db):
    for tr in genes[query.chrom].find(query.start, query.end):
        try:
            ensid = trans_db[tr.value['geneid']]
            print "%s\t%s" % (tr.value['geneid'], ensid)
        except KeyError:
            print >> sys.stderr, "%s\t%s" % (tr.value['geneid'], "unknown")

def main():
    if len(sys.argv) < 5:
        print >> sys.stderr, \
            'Usage: get_ens_name.py bed_file result_file coord_file trans_db'
        raise SystemExit

    bedfile = sys.argv[1]
    result_file = sys.argv[2]
    coord_file = sys.argv[3]
    trans_db = cPickle.load(open(sys.argv[4]))

    genes = {}
    exon_bins = {}

    print >> sys.stderr, 'Reading', bedfile

    for n, exons in enumerate(parse_models(bedfile), start=1):
        add_intervals(genes, exons)

    print >> sys.stderr, 'Reading', coord_file
    bins = csv.reader(open(coord_file), dialect='excel-tab')
    for row in bins:
        bin_name = row[0]+"_"+row[1]
        exon_bins[bin_name] = ExonQuery(row[2], int(row[3]), int(row[4]))

    print >> sys.stderr, 'Reading', result_file
    results = csv.reader(open(result_file), dialect='excel-tab')
    for row in results:
        agg_name = row[0]+"_"+row[1]
        query = exon_bins[agg_name]
        find_id(genes, query, trans_db)

if __name__=='__main__':
    main()
