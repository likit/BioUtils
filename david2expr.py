'''The script reads in a list of genes from an output from DAVID and reads in
expression data from a text file containing corresponding genes and report
an expression level of each gene in DAVID list to standard output.

Both gene expression and DAVID files should be in a comma-delimited format.

'''

import sys

david_file = sys.argv[1]
expr_file = sys.argv[2]

genes = {}

print >> sys.stderr, 'Reading %s...', david_file

for line in open(expr_file):
    cols = line.strip().split(',')
    geneid = cols[0]
    exprs = '\t'.join(cols[1:5])
    if geneid not in genes:
        genes[geneid] = exprs
    else:
        raise KeyError('duplicated gene ID')

print >> sys.stderr, 'Reading %s...', expr_file

for line in open(david_file):
    geneid = line.strip().split(',')[0]
    if geneid in genes:
        print '%s\t%s' % (geneid, genes[geneid])
