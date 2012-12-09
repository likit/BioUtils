import sys

import networkx as nx
from pygr import seqdb, cnestedlist, sequtil
import matplotlib.pyplot as plt

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

        yield AlignObj(target, query, t_start, t_end,
                                        q_start, q_end)

def add_alignment(aligndb, al, targets, queries):
    if al.t_end < al.t_start:
        al.t_start, al.t_end = al.t_end, al.t_start
        t_ival = -targets[al.target][al.t_start:al.t_end]
    else:
        t_ival = targets[al.target][al.t_start:al.t_end]
    if al.q_end < al.q_start:
        al.q_start, al.q_end = al.q_end, al.q_start
        q_ival = -queries[al.query][al.q_start:al.q_end]
    else:
        q_ival = queries[al.query][al.q_start:al.q_end]

    aligndb[t_ival] += q_ival


def write_sequence(source, graph, targets, queries, ofile):
    visited_nodes = set()
    max_length = 0
    for node in nx.algorithms.dfs_preorder_nodes(graph, source):
        seq = targets[node] if node in targets else queries[node]

        if len(seq) > max_length: max_length = len(seq)

        sequtil.write_fasta(ofile, str(seq), 60, node)
        visited_nodes.add(node)

    return visited_nodes, max_length

def main():
    if len(sys.argv) < 4: raise SystemExit

    print >> sys.stderr, 'Reading sequence databases...'
    queries = seqdb.SequenceFileDB(sys.argv[1])
    targets = seqdb.SequenceFileDB(sys.argv[2])
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

    for al in parse_alignments(align_file):
        aligndb += targets[al.target]
        add_alignment(aligndb, al, targets, queries)

    print >> sys.stderr, 'Building the alignment database...'
    aligndb.build()

    print >> sys.stderr, 'Constructing alignment graphs...'
    graph = nx.Graph()
    for target in targets:
        try:
            sub_ival = targets[target]
            for src, dest, edge in aligndb[sub_ival].edges():
                source = repr(src).split('[')[0].lstrip('-')
                destination = repr(dest).split('[')[0].lstrip('-')
                graph.add_edge(source, destination)
        except KeyError:
            pass

    # nx.draw(graph)
    # plt.show()
    # print graph.nodes()
    logfile = open('assemgraph.log', 'w')
    visited_nodes = set()
    cluster_no = 0
    for node in graph.nodes():
        if node not in visited_nodes:
            filename = 'cluster_%d' % cluster_no
            ofile = open(filename, 'w')
            print >> sys.stderr, 'Writing %s to a file...' % filename,
            vnodes, max_length = (write_sequence(node, graph, targets,
                                                        queries, ofile))
            visited_nodes.update(vnodes)
            for n in vnodes:
                size = len(targets[n]) if n in targets else len(queries[n])
                print >> logfile, '%s\t%s\t%d' % (filename, n, size)
            ofile.close()
            print >> sys.stderr, '\ttotal nodes = %d' % len(vnodes)

            cluster_no += 1

    logfile.close()


if __name__=='__main__':
    main()
