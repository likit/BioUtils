import sys

import networkx as nx
from pygr import seqdb, sequtil

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
        score = float(cols[11])

        if score < MIN_SCORE: continue

        yield query, target

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
    if len(sys.argv) < 4:
        print >> sys.stderr, \
                'Usage: assemgraph2.py <query file> <target file>' + \
                ' <blast9 file> [min score=200]'
        raise SystemExit
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

    print >> sys.stderr, 'Constructing alignment graphs...'
    graph = nx.Graph()
    for c, (query, target) in enumerate(parse_alignments(align_file)):
        graph.add_edge(query, target)

        if c % 100 == 0:
            print >> sys.stderr, '...', c

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
