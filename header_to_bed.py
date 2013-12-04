'''Creates BED file from SAM header.'''

import sys

if __name__=='__main__':
    try:
        infile = open(sys.argv[1])
    except:
        infile = sys.stdin
    for line in infile:
        if line.startswith('@SQ'):
            seqname, seqlen = line.split()[1:]
            seqname = seqname.split(':')[1]
            seqlen = int(seqlen.split(':')[1])
            print '%s\t0\t%d' % (seqname, seqlen)

