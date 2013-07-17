'''Split interleaved paired-end reads to two files.'''

import sys

if __name__=='__main__':
    op1 = open('reads1.fq', 'w')
    op2 = open('reads2.fq', 'w')
    fp = open(sys.argv[1])

    n = 0
    while True:
        try:
            line = fp.readline()
            if line == "":
                break
            n += 1
        except:
            break
        else:
            if line.endswith('/1\n'):
                print >> op1, line,
                print >> op1, fp.readline(),
                print >> op1, fp.readline(),
                print >> op1, fp.readline().strip()
            if line.endswith('/2\n'):
                print >> op2, line,
                print >> op2, fp.readline(),
                print >> op2, fp.readline(),
                print >> op2, fp.readline().strip()
            if (n % 100000) == 0:
                print >> sys.stderr, '... %d' % n
    op1.close()
    op2.close()

    if n == 0:
        print >> sys.stderr, "Cannot open file %s, or the file is empty." \
                % sys.argv[1]
