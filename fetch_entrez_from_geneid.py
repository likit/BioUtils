'''Selects protein sequences from NCBI that are in a list
from Geisha text file.

Output is written to standard output.

'''

import os
import sys
import time
from Bio import SeqIO, Entrez

def parse(infile):
    '''Return a set of gene IDs from an input file.'''

    for line in open(infile):
        geneid = line.split()[0]

        yield geneid

def fetch(geneid):
    print >> sys.stderr, 'fetching.. gene ID: %s' % geneid
    handle = Entrez.efetch(db='gene', retmode='xml', id=geneid)
    xmldata = Entrez.read(handle)
    product = xmldata[0]['Entrezgene_locus'][0]\
                ['Gene-commentary_products'][0]
    prodtype = product['Gene-commentary_type'].attributes['value']
    print >> sys.stderr, 'product type = %s' % (prodtype)

    seq_gi = xmldata[0]['Entrezgene_locus'][0]\
                ['Gene-commentary_products'][0]\
                ['Gene-commentary_seqs'][0]\
                ['Seq-loc_whole']['Seq-id']\
                ['Seq-id_gi']

    handle = Entrez.efetch(db='nucleotide', retmode='text',
                            rettype='fasta', id=seq_gi)

    seq = SeqIO.read(handle, 'fasta')
    return seq

def main():
    infile = sys.argv[1]
    Entrez.email = sys.argv[2]

    outfile = os.path.splitext(infile)[0] + ".fa"
    records = []

    for geneid in parse(infile):
        try:
            records.append(fetch(geneid))
        except:
            print >> sys.stderr, 'Cannot retrieve a sequence'
            continue

        time.sleep(3)

    SeqIO.write(records, outfile, 'fasta')

    print >> sys.stderr, 'Total sequences = %d' % len(records)

if __name__=='__main__':
    main()
