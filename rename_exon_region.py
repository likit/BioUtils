'''The script rename exon regions from flattened GFF file
created by dexseq_prepare_annotation.py to a shorter name.
The input is an output from dexseq_count.py 
A new name is quite arbitrary, depending on the order of
a region in the input file and a given prefix.

Author: Likit Preeyanon
Email: preeyano@msu.edu

'''

import sys

def rename_count(afile, ref_file):
    refs = {}
    for line in open(ref_file):
        region, name = line.strip().split('\t')
        refs[name] = region

    for line in open(afile):
        if line.startswith('_'):
            continue
        name, value = line.strip().split()
        name, exon_id = name.split(':')
        print '%s:%s\t%s' % (refs[name], exon_id, value)


def rename_gff(afile, prefix):
    region_id = 1
    for line in open(afile):
        cols = line.strip().split('\t')
        if cols[2] == 'aggregate_gene':
            region = prefix + '-' + str(region_id)
            ori_geneid = cols[-1].split()[-1].strip('"')
            print >> sys.stderr, '%s\t%s' % (region, ori_geneid)

            cols[-1] = 'gene_id "%s"' % (region)
            region_id += 1
        else:
            transcript, exon_part, geneid = cols[-1].split('; ')
            cols[-1] = '%s; %s; gene_id "%s"' % (transcript, exon_part, region)
        print '\t'.join(cols)


if __name__=='__main__':
    try:
        input_file_type = sys.argv[1]
        afile = sys.argv[2]
    except IndexError:
        print >> sys.stderr, 'Need more arguments.'
        print >> sys.stderr, 'Usage: rename_exon_region.py <type GFF/count> ' +\
                                '<GFF/count file> <prefix>'
    else:
        if input_file_type == 'count':
            ref_file = sys.argv[3]
            rename_count(afile, ref_file)
        elif input_file_type == 'gff':
            prefix = sys.argv[3]
            rename_gff(afile, prefix)
        else:
            print >> sys.stderr, 'Unrecognized file type.'
