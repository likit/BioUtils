'''Reads a BED file and a WIG file and report the highest coverage of
an interval specified in BED file.
An interval can be exon, gene, peak and etc.
WIG file and BED file should store data from the same chromosome.
Output is written to standard output in tab-delimited format.

'''
from bx.intervals.intersection import Interval, IntervalTree
import sys

USAGE = '''python peak-coverage-csv.py <BED file> <WIG file>'''

def parse_wig(wigfile):
    '''Parse positions and coverages from WIG file and return
    an interval with maximum coverage
    
    '''

    coverages = set()
    start = None
    end = None

    print >> sys.stderr, 'Parsing %s' % wigfile

    for line in open(wigfile):
        try:
            pos, cov = line.strip().split()
            pos = int(pos)
            cov = float(cov)
        except ValueError:
            continue

        if not start: # begining of the file
            start = pos
            end = start
            coverages.add(cov)
        else:
            if pos != end + 1: # positions not connected

                yield Interval(start, end, value={'max_cov':max(coverages)})

                start = pos
                end = start
                coverages = set([cov]) # reset the coverages
            else:
                end += 1 # move to the next position
                coverages.add(cov)

    yield Interval(start, end, value={'max_cov':max(coverages)})

def parse_bed(bedfile):
    '''Return start and end position of an interval
    An interval can be exon, gene, peak etc.
    '''

    print >> sys.stderr, 'Parsing %s' % bedfile

    for line in open(bedfile):
        try:
            chrom, start, end = line.strip().split()[:3]
            start, end = int(start), int(end)
        except:
            continue
        else:
            yield chrom, start, end

def get_coverage(intersects):
    '''Return average coverage of multiple intersections.
    Return None if there are no intersections.

    '''

    cov = 0
    if not intersects: # no intersections
        return None

    for intr in intersects:
        cov += intr.value['max_cov']

    if len(intersects) > 1:
        print >> sys.stderr, len(intersects)
    return cov/len(intersects)

def report(intersecter, bedfile):
    '''Write results to standard output in tab delimited format as
    1) Chromosome
    2) Start
    3) End
    4) Coverage

    '''

    for chrom, start, end in parse_bed(bedfile):
        '''To cover a coverage of a single position,
        start and end positions are extended one position each.

        '''
        start = start - 1 if start - 1 > 0 else start
        end += 1

        avg_cov = get_coverage(intersecter.find(start, end))
        if avg_cov:
            print '%s\t%d\t%d\t%.5f' % (chrom, start, end, avg_cov)

def main(argv):
    bedfile = argv[1]
    wigfile = argv[2]
    intersecter = IntervalTree()
    for peak in parse_wig(wigfile):
        intersecter.insert_interval(peak)

    report(intersecter, bedfile)

if __name__=='__main__':
    try:
        main(sys.argv)
    except:
        print >> sys.stderr, USAGE
