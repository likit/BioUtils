import sys
from collections import namedtuple

Read = namedtuple('Read', ['name','qual','seq'])

read1 = None
fstmate = open(sys.argv[1], 'w')
scndmate = open(sys.argv[2], 'w')
unpaired = open(sys.argv[3], 'w')
for line in sys.stdin:
    items = line.split('\t')
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
