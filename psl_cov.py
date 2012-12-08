'''Select alignments with >= a given coverage.'''

import sys
MIN_COV = 0.6
fails = 0

for n, line in enumerate(open(sys.argv[1]), start=1):
    cols = line.strip().split()
    qSize = int(cols[10])
    qStart = int(cols[11])
    qEnd = int(cols[12])
    coverage = float(qEnd - qStart)/qSize
    if coverage >= MIN_COV:
        print line.strip()
    else:
        fails += 1

print >> sys.stderr, 'Total alignment = %d\nTotal failed alignments = %d' % \
                        (n, fails)
