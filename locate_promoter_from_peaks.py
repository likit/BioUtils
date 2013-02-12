'''The script reads data download from UCSC genome browser
and a peak file containing a gene name, TSS, start and end and
determine is a given peak located in promoter region of a particular
gene.

'''

import sys

class Gene(object):
    def __init__(self, chrom, strand, tss_start, tss_end, refname):
        self.chrom = chrom
        self.strand = strand
        self.tss_start = tss_start
        self.tss_end = tss_end
        self.refname = refname

def build_genes_db(dbfile):
    genes_db = {}
    for line in open(dbfile):
        cols = line.split()
        if len(cols) < 6:
            print line
            sys.exit()

        refname = cols[0]
        chrom = cols[1]
        strand = cols[2]
        tss_start = int(cols[3])
        tss_end = int(cols[4])
        name = cols[5]
        gene = Gene(chrom, strand, tss_start, tss_end, refname)
        if name not in genes_db:
            genes_db[name] = [gene]
        else:
            genes_db[name].append(gene)

    return genes_db

def is_promoter(infile, genes_db):
    print "#Gene_name\trefseq_id\tchrom\tstrand\tOld_TSS_start\t" + \
            "TSS_start\tTSS_end\tpeak_start\tpeak_end\t" + \
            "distance\tpromoter\tremark"

    for line in open(infile):
        try:
            name, tss, start, end = line.split()
        except ValueError:
            print line
            sys.exit()

        start = int(start)
        end = int(end)
        try:
            gene = genes_db[name]
        except KeyError:
            continue

        for isoform in gene:
            if isoform.strand == '+':
                dist = isoform.tss_start - start
            else:
                dist = start - isoform.tss_end
            promoter = "YES" if (dist > 0 and dist <= 10000) else "NO"
            tss = int(tss)
            remark = "*" if tss != isoform.tss_start else ""

            print "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s" % (
                    name, isoform.refname, isoform.chrom, isoform.strand,
                    tss, isoform.tss_start, isoform.tss_end, start, end,
                    dist, promoter, remark)

def main():
    dbfile = sys.argv[1]
    infile = sys.argv[2]
    genes_db = build_genes_db(dbfile)
    is_promoter(infile, genes_db)

if __name__=='__main__':
    main()
