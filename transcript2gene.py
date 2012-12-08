'''The script reads gene models in BED format and creates BED file with only
gene name and its coordinate that covers all isoforms to be used with BEDTools
for DGE analysis.

'''

import sys
import csv

genes = {}
gene_list = []

def parse_bed(bed_file):
    '''Reads alignments from BED format and creates
    exon objects from a transcript.

    '''
    reader = csv.reader(open(bed_file), dialect='excel-tab')
    for row in reader:
        chrom = row[0]
        chrom_start = int(row[1])
        chrom_end = int(row[2])
        geneid = row[3].split('.')[0]
        if geneid not in genes:
            gene_list.append(geneid)
            genes[geneid] = [chrom_start, chrom_end, chrom]
        else:
            if chrom_start < genes[geneid][0]:
                genes[geneid][0] = chrom_start
            if chrom_end < genes[geneid][1]:
                genes[geneid][1] = chrom_end

def main():
    bed_file = sys.argv[1]
    parse_bed(bed_file)
    for gene in gene_list:
        start, end, chrom = genes[gene]
        print '%s\t%d\t%d\t%s' % (chrom, start, end, gene)

if __name__=='__main__':
    main()
