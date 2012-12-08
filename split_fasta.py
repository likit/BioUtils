'''Split sequences in FASTA format based on a given keyword
in a sequence title
'''

import sys
from Bio import SeqIO

keyword = sys.argv[2]

output_file = sys.argv[3]

select = []

for n, sequence in enumerate(SeqIO.parse(sys.argv[1], "fasta")):
    if keyword in sequence.id:
        select.append(sequence)
    if n%1000==0:
        print >> sys.stderr, '...', n

SeqIO.write(select, output_file, "fasta")
