import sys

import networkx as nx
from pygr import seqdb, cnestedlist, sequtil
import matplotlib.pyplot as plt

MIN_SCORE = 200

class AlignObj(object):
    def __init__(self, target, query, t_start, t_end,
                                        q_start, q_end):
        self.target = target
        self.query = query
        self.t_start = t_start
        self.t_end = t_end
        self.q_start = q_start
        self.q_end = q_end


def parse_alignments(align_file):
    for line in align_file:
        if line.startswith('#'): continue
        cols = line.strip().split()
        query = cols[0]
        target = cols[1]
        q_start = int(cols[6])
        q_end = int(cols[7])
        t_start = int(cols[8])
        t_end = int(cols[9])
        score = float(cols[11])

        if score < MIN_SCORE: continue
        yield AlignObj(target, query, t_start, t_end,
                                        q_start, q_end)

def add_alignment(aligndb, al, targets, queries):
    # try:
    if al.t_end < al.t_start:
        al.t_start, al.t_end = al.t_end, al.t_start
        t_ival = -targets[al.target][al.t_start - 1:al.t_end]
    else:
        t_ival = targets[al.target][al.t_start - 1:al.t_end]
    if al.q_end < al.q_start:
        al.q_start, al.q_end = al.q_end, al.q_start
        q_ival = -queries[al.query][al.q_start - 1:al.q_end]
    else:
        q_ival = queries[al.query][al.q_start - 1:al.q_end]
    # except IndexError:
    #     print >> sys.stderr, al.target, al.query, al.t_start, al.t_end, al.q_start, al.q_end
    #     raise IndexError

    aligndb[t_ival] += q_ival


def write_sequence(source, graph, targets, queries, ofile1, ofile2):
    visited_nodes = set()
    max_length = 0
    for node in nx.algorithms.dfs_preorder_nodes(graph, source):
        seq = targets[node] if node in targets else queries[node]

        if len(seq) > max_length: max_length = len(seq)

        if node in targets:
            sequtil.write_fasta(ofile1, str(seq), 60, node)
        else:
            sequtil.write_fasta(ofile2, str(seq), 60, node)

        visited_nodes.add(node)

    return visited_nodes, max_length

def main():
    if len(sys.argv) < 4: raise SystemExit
    try:
        MIN_SCORE = float(sys.argv[5])
    except IndexError:
        pass

    print >> sys.stderr, 'Reading sequence databases...'
    queries = seqdb.SequenceFileDB(sys.argv[1])
    targets = seqdb.SequenceFileDB(sys.argv[2])
    print >> sys.stderr, len(queries), len(targets)
    try:
        align_file = open(sys.argv[3])
    except IOError as e:
        print >> sys.stderr, 'Error: check alignment file.'
        raise e

    aligndb = cnestedlist.NLMSA('alignment', mode='memory',
                                                pairwiseMode=True)

    print >> sys.stderr, 'Adding sequences to an alignment database...'
    # for n, target in enumerate(targets):
    #     aligndb += targets[target]
    #     if n % 1000 == 0: print >> sys.stderr, '...', n

    target_list = set()

    for c, al in enumerate(parse_alignments(align_file)):
        aligndb += targets[al.target]
        target_list.add(al.target)
        add_alignment(aligndb, al, targets, queries)
        if c % 100 == 0: print >> sys.stderr, '...', c

    print >> sys.stderr, 'Building the alignment database...'
    aligndb.build()

    print >> sys.stderr, 'Constructing alignment graphs...'
    graph = nx.Graph()
    for c, target in enumerate(target_list):
        try:
            sub_ival = targets[target]
            for src, dest, edge in aligndb[sub_ival].edges():
                source = repr(src).split('[')[0].lstrip('-')
                destination = repr(dest).split('[')[0].lstrip('-')
                graph.add_edge(source, destination)
        except KeyError:
            pass
        if c % 100 == 0: print >> sys.stderr, '...', c

    # nx.draw(graph)
    # plt.show()
    # print graph.nodes()
    logfile = open('assemgraph.log', 'w')
    visited_nodes = set()
    cluster_no = 0
    for node in graph.nodes():
        if node not in visited_nodes:
            filename1 = 'cluster_%d_targets' % cluster_no
            filename2 = 'cluster_%d_queries' % cluster_no
            ofile1 = open(filename1, 'w')
            ofile2 = open(filename2, 'w')
            print >> sys.stderr, \
                    'Writing cluster %d to a file...' % cluster_no,
            vnodes, max_length = (write_sequence(node, graph, targets,
                                                queries, ofile1, ofile2))
            visited_nodes.update(vnodes)
            for n in vnodes:
                size = len(targets[n]) if n in targets else len(queries[n])
                print >> logfile, 'cluster_%d\t%s\t%d' % (cluster_no, n, size)
            ofile1.close()
            ofile2.close()
            print >> sys.stderr, '\ttotal nodes = %d' % len(vnodes)

            cluster_no += 1

    print >> logfile, '***finished***'
    logfile.close()


if __name__=='__main__':
    main()
