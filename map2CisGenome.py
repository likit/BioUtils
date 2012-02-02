'''The script converts Bowtie format to CisGenome format.

Author: Likit Preeyano, preeyano@msu.edu.

'''

import sys, csv

def convert(filename):
    reader = csv.reader(open(filename), dialect='excel-tab')
    writer = csv.writer(sys.stdout, dialect='excel-tab')
    for row in reader:
        strand = 'F' if row[1]=='+' else 'R'
        chrom = row[2]
        start = row[3]
        writer.writerow([chrom, start, strand])


if __name__=='__main__':
    convert(sys.argv[1])
