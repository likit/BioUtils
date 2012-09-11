'''Reads a BED file and a WIG file and report the highest coverage of
each location specified in BED file in CSV format (tab-delimited).

The script is intended to use in Chip-Seq analysis.

Output is written to standard output.

'''
import bx
import sys
from collections import namedtuple

USAGE = '''python peak-coverage-csv.py <BED file> <WIG file>'''

Peak = namedtuple('Peak', ['start', 'end', 'max_cov'])

def parse_wig(wigfile):
    coverages = set()
    start = None
    end = None
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

                yield Peak(start, end, max(coverages))

                start = pos
                end = start
                coverages = set([cov]) # reset the coverages
            else:
                end += 1 # move to the next position
                coverages.add(cov)

if __name__=='__main__':
    for peak in parse_wig(sys.argv[1]):
        print peak.start, peak.end, peak.max_cov
