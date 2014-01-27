'''Extracts paired and single-end reads from BAM file.

BAM file has to be sorted by name. The following commands
can be used to sort and extract reads.

samtools sort -n unmapped.bam unmapped.sorted
samtools view unmapped.sorted.bam | python sam2fq.py r1.fq r2.fq un.fq

where r1.fq, r2.fq and un.fq are output files for
first mate, second mate and unpaired reads respectively.
'''

import sys
from collections import namedtuple

Read = namedtuple('Read', ['name','qual','seq'])

read1 = None
fstmate = open(sys.argv[1], 'w')
scndmate = open(sys.argv[2], 'w')
unpaired = open(sys.argv[3], 'w')
for line in sys.stdin:
    items = line.strip().split('\t')
    name, qual, seq = items[0], items[10], items[9]
    if not read1:
        read1 = Read(name, qual, seq)
        continue
    else:
        read2 = Read(name, qual, seq)

    if read1.name == read2.name:
        print >> fstmate, '@%s\n%s\n+\n%s' % (read1.name, read1.seq, read1.qual)
        print >> scndmate, '@%s\n%s\n+\n%s' % (read2.name, read2.seq, read2.qual)
        read1 = None
    else:
        print >> unpaired, '@%s\n%s\n+\n%s' % (read1.name, read1.seq, read1.qual)
        read1 = read2
        read2 = None

if read1:
    print >> unpaired, '@%s\n%s\n+\n%s' % (read1.name, read1.seq, read1.qual)
