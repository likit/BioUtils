'''Reads sequence IDs or gene names from a text file and retrieve sequences
from Entrez Db.

'''

import sys
from Bio import Entrez, SeqIO

def fetch(queries, email, db="protein",
            rettype="fasta", retmode="text"):
    for qid in queries:
        handle = Entrez.efetch(db=db, rettype=rettype,
                                retmode=retmode, id=qid)
        seq_record = SeqIO.read(handle, rettype)
        SeqIO.write(seq_record, sys.stdout, 'fasta')
        print >> sys.stderr, '%s...' % seq_record.description[:45]
        handle.close()

def read_id(id_file):
    queries = []
    for line in open(sys.argv[1]):
        queries.append(line.strip())

    return queries

if __name__=='__main__':
    try:
        id_file = sys.argv[1]
        Entrez.email = sys.argv[2]
    except IndexError:
        print >> sys.stderr, 'Usage: python fetch_seq.py <id file> <email>'
        raise SystemExit
    else:
        queries = read_id(id_file)
        fetch(queries, Entrez.email)
