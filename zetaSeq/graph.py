# !/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri May 27 09:29:44 2022

He who loves to comment his code is unlikely to have bad luck.
@author: Song, Zewei
@contact: songzewei@genomics.cn
"""

# %% Oh all the functions

import random
from zetaSeq import io as seqIO
import json
import networkx as nx
from collections import OrderedDict


# Generate a random sequence with length = l
def random_sequence(l):
    n = ('A', 'T', 'C', 'G')
    s = ''
    for i in range(l):
        s += random.sample(n, 1)[0]
    return s


# Reverse compliment a sequence
def rc(seq):
    d = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    seq_rc = ''.join([d[i] for i in seq])[::-1]
    return seq_rc


# Generate a random fragment with a given length at (minlen, maxlen)
def random_fragment(seq, minlen, maxlen):
    start = int(random.random() * len(seq))  # Pick a random start
    r = random.sample(range(minlen, maxlen, 1), 1)[0]  # Pick a random length
    # strand = random.sample((True,False), 1)[0]
    strand = True
    if start + r > len(seq):
        if strand:
            return seq[start:]
        else:
            return rc(seq[start:])
    else:
        if strand:
            return seq[start:start + r + 1]
        else:
            return rc(seq[start:start + r + 1])


# Read a PAF alignment file into dictionary
def load_paf(input_paf):
    alignment = []
    with open(input_paf, 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            line[1] = int(line[1])
            line[2] = int(line[2])
            line[3] = int(line[3])
            line[6] = int(line[6])
            line[7] = int(line[7])
            line[8] = int(line[8])
            line[9] = int(line[9])
            line[10] = int(line[10])
            if line[0] != line[5]:
                alignment.append(
                    {'query': line[0], 'qlen': line[1], 'qstart': line[2], 'qend': line[3], 'strand': line[4],
                     'target': line[5], 'tlen': line[6], 'tstart': line[7], 'tend': line[8],
                     'match': line[9], 'align': line[10]})
    return alignment


# Return a type key for overlap():
def type_key(item, t=300):
    s = [('qs', item['qstart']), ('ts', item['tstart']), \
         ('qle', item['qlen'] - item['qend']), ('tle', item['tlen'] - item['tend']), \
         ('t', t), ('2t', t * 2)]
    s.sort(key=lambda i: i[1])
    a = sorted([s[0], s[1]], key=lambda x: x[0])
    b = sorted([s[4], s[5]], key=lambda x: x[0])
    s[0] = a[0]
    s[1] = a[1]
    s[4] = b[0]
    s[5] = b[1]
    return tuple([x[0] for x in s])


# Return overlap pairs from input alignments
# Asign one PAF line to an alignment type (TYPE 3 or TYPE 4), with a tail threshold (300 at default)
# Return ((query, strd), (target, strd), match, (query_end, target_start), (target'_end, query'_start))
def overlap(item, t):
    overlap_type = {('qle', 'ts', 't', '2t', 'qs', 'tle'): ('+', 3),
                    ('qs', 'tle', 't', '2t', 'qle', 'ts'): ('+', 4),
                    ('qle', 'tle', 't', '2t', 'qs', 'ts'): ('-', 3),
                    ('qs', 'ts', 't', '2t', 'qle', 'tle'): ('-', 4)}
    typ = overlap_type.get(type_key(item, t), '')
    qes = item['qend'] - item['qstart']
    tes = item['tend'] - item['tstart']
    if typ == ('+', 3):
        return ((item['query'], '+'), (item['target'], '+'), item['match'],
                (qes + item['qlen'] - item['qend'], item['tstart']),
                (tes + item['tstart'], item['qlen'] - item['qend']))
    elif typ == ('+', 4):
        return ((item['target'], '+'), (item['query'], '+'), item['match'],
                (tes + item['tlen'] - item['tend'], item['qstart']),
                (qes + item['qstart'], item['tlen'] - item['tend']))
    elif typ == ('-', 3):
        return ((item['query'], '+'), (item['target'], '-'), item['match'],
                (qes + item['qlen'] - item['qend'], item['tlen'] - item['tend']),
                (tes + item['tlen'] - item['tend'], item['qlen'] - item['qend']))
    elif typ == ('-', 4):
        return ((item['target'], '-'), (item['query'], '+'), item['match'],
                (item['tend'], item['qstart']),
                (item['qend'], item['tstart']))
    else:
        return None


# Return a list contains all overlaps
def overlap_list(paf_file, t=300, id=0.99, m=1000):
    alignment = load_paf(paf_file)
    ovlp = [overlap(item, t) for item in alignment
            if item['match'] / item['align'] > id
            and item['align'] >= m]
    ovlp = [i for i in ovlp if i != None]  # Keep only eligible overlap patterns
    ovlp = set(ovlp)
    return list(ovlp)


# Return the reverse complement of the ovrlap pair
def rc_overlap(item):
    strd = {'+': '-', '-': '+'}
    rc_ovlp = ((item[1][0], strd[item[1][1]]), (item[0][0], strd[item[0][1]]), item[2],
               item[4], item[3])
    return rc_ovlp


# Return the reverse complement of the given node
def rc_node(item):
    strd = {'+': '-', '-': '+'}
    return (item[0], strd[item[1]])


# Given a node n in g, return all nodes connected to it (thus the subgraph of n)
def max_subgraph(g, n):
    nodes = set([i for i in nx.all_neighbors(g, n)] + [n])
    added_nodes = len(nodes)
    while added_nodes > 0:
        l = [len(nodes), 0]
        new_nodes = []
        for item in nodes:
            new_nodes += [i for i in nx.all_neighbors(g, item)]
        nodes = nodes.union(set(new_nodes))
        l[1] = len(nodes)
        added_nodes = l[1] - l[0]
    return nodes


# Given a graph g, return a list of its max_subgraph
def subgraph(g):
    nodes = set([i for i in g.nodes])  # get all nodes in the graph g
    subg = []  # list for all subgraphs
    while nodes:  # while there is nodes unassigned
        new_nodes = max_subgraph(g, list(nodes)[0])
        subg.append(nx.subgraph(g, new_nodes))
        nodes = nodes.difference(new_nodes)
    return subg


# Given a graph, remove all vertex with indegree > 1 or outdegree > 1
# aka remove all branching edges.
def remove_branch_edges(g):
    for node in g.nodes:
        successors = [i for i in g.successors(node)]
        predecessors = [i for i in g.predecessors(node)]
        if len(successors) > 1:
            for item in successors:
                g.remove_edge(node, item)
        if len(predecessors) > 1:
            for item in predecessors:
                g.remove_edge(item, node)
    return g


# Transitive reduction of the given graph
# Now that we omit all directed cyclic graph, will add later
def reduce_graph(g):
    # remove triad 030T (three nodes)
    if nx.is_directed_acyclic_graph(g):
        trg = nx.transitive_reduction(g)
        # add node and edge values to the new graph
        trg.add_nodes_from(g.nodes(data=True))
        trg.add_edges_from((u, v, g.edges[u, v]) for u, v in trg.edges)
        return trg
    else:
        return g


# Return a sequence with Given label and strand
def get_sequence(node, sequence_dict):
    if node[1] == '+':
        return sequence_dict[node]
    elif node[1] == '-':
        return rc(sequence_dict[(node[0], '+')])
    else:
        return None


# Return a ovlp_dict object for next_sequence()
def overlap_dict(ovlp):
    ovlp_dict = {(i[0], i[1]): i for i in ovlp}
    for item in ovlp:
        rc_item = rc_overlap(item)
        ovlp_dict[(rc_item[0], rc_item[1])] = rc_item
    return ovlp_dict


# Return a sequence_dict object for next_sequence()
def sequence_dict(seq_file):
    c = {(i[0], '+'): i[1] for i in seqIO.sequence(seq_file)}
    return c


# Reture an overlapped sequences from previous sequence and current edge
def next_sequence(prev_seq, tail_node, next_node, ovlp_dict, sequence_dict):
    ovlp = ovlp_dict[(tail_node, next_node)]  # get the next edge
    cutoff_tail = len(prev_seq) - ovlp[3][0]
    cutoff_head = ovlp[3][1]
    next_seq = get_sequence(next_node, sequence_dict)
    seq = prev_seq[:cutoff_tail] + next_seq[cutoff_head:]
    return seq


# Return the overlap sequence from a garph with chained nodes
def ovlp_sequence(g, ovlp_dict, seq_dict):
    first_node = [i for i in g.nodes if g.in_degree(i) == 0]
    current_node = first_node[0]
    label = current_node[0] + '[' + current_node[1] + ']'
    current_seq = get_sequence(current_node, seq_dict)
    while list(g.successors(current_node)):  # Continue if this node has a successor node
        next_node = list(g.successors(current_node))[0]
        label += ':' + next_node[0] + '[' + next_node[1] + ']'
        current_seq = next_sequence(current_seq, current_node, next_node, ovlp_dict, seq_dict)
        current_node = next_node
    return [label, current_seq]


# Build an overlap layouot graph with a given overlap_list from overlap_list()
def overlap_layout_graph(ovlp, log=False):
    # Prepare for searching dictionary
    query_ovlp = {}  # ovlp dictionary using query as keys
    target_ovlp = {}  # ovlp dictionary using target as keys
    rc_ovlp = []  # list for rc of all ovlp
    if log: print('Found {0} edges'.format(len(ovlp)))
    for item in ovlp:
        query = item[0]  # I'm query
        target = item[1]  # I'm target
        query_ovlp[query] = query_ovlp.get(query, []) + [item]
        target_ovlp[target] = target_ovlp.get(target, []) + [item]
        rc_item = rc_overlap(item)  # I'm the reverse compliment of myself
        rc_ovlp.append(rc_item)
        query = rc_item[0]
        target = rc_item[1]
        query_ovlp[query] = query_ovlp.get(query, []) + [rc_item]
        target_ovlp[target] = target_ovlp.get(target, []) + [rc_item]
    query_ovlp = {i: set(j) for i, j in query_ovlp.items()}  # set all values
    target_ovlp = {i: set(j) for i, j in target_ovlp.items()}  # set all values
    ovlp = ovlp + rc_ovlp  # all possible ovlp considering rc cases
    ovlp.sort(key=lambda i: i[2], reverse=True)  # reverse sort by matching length
    ovlp = OrderedDict.fromkeys(ovlp)  # set using OrderedDict

    # Add ovlp to DiGraph(), creating new graph is necessary
    g_list = []  # List for overlap graphs
    if not g_list:  # if g_list is empty, add the first ovlp to the graph
        g_list.append(nx.DiGraph())  # add the first DiGraph() object
        adding_ovlp = list(ovlp.keys())[0]  # add the first overlap pair
        if log: print('The first to add is {0}'.format(adding_ovlp))
        g_list[-1].add_edge(adding_ovlp[0], adding_ovlp[1], match=adding_ovlp[2])
        ovlp.pop(adding_ovlp)  # remove the added overlap pair
        ovlp.pop(rc_overlap(adding_ovlp))  # remove its rc too
        last_added_nodes = set([adding_ovlp[0], adding_ovlp[1]])
        if log: print('Current length of ovlp: {0}'.format(len(ovlp)))

    # Add overlap pair from ovlp to current graphs, until ovlp is empty
    while ovlp:
        added_nodes = set([])  # set added_nodes as a blank set
        # Loop for nodes in the latest graph, add edge start from these nodes if no conflict exist
        for node in last_added_nodes:
            if node in query_ovlp:  # check for nodes in the query position
                for edge in query_ovlp[node]:
                    if edge in ovlp:
                        # check if target exist in current graph in rc form:
                        if not g_list[-1].has_node(rc_node(edge[1])):
                            adding_ovlp = edge
                            g_list[-1].add_edge(edge[0], edge[1], match=edge[2])
                            ovlp.pop(adding_ovlp)  # remove added edge fron ovlp
                            ovlp.pop(rc_overlap(adding_ovlp))  # remove its rc too
                            added_nodes.add(edge[1])  # add the new node to the added list
            if node in target_ovlp:  # check for nodes in the target position
                for edge in target_ovlp[node]:
                    if edge in ovlp:
                        # check if target exist in current graph in rc form:
                        if not g_list[-1].has_node(rc_node(edge[0])):
                            adding_ovlp = edge
                            g_list[-1].add_edge(edge[0], edge[1], match=edge[2])
                            ovlp.pop(adding_ovlp)  # remove added edge from ovlp
                            ovlp.pop(rc_overlap(adding_ovlp))  # rmove its rc too
                            added_nodes.add(edge[0])  # add the new node to the added list
        last_added_nodes = added_nodes  # assign added_nodes back to last_added_nodes
        if log: print('Added nodes this round: {0}'.format(added_nodes))
        if log: print('{0} graphs so far'.format(len(g_list)))

        # Add a new graph if len(added_nodes) == 0, aka no new node can be added
        if len(added_nodes) == 0:  # if no new nodes can be added to current graph
            if log: print('Added a new graph')
            g_list.append(nx.DiGraph())  # add a new DiGraph to the list
            adding_ovlp = list(ovlp.keys())[0]  # add the first overlap pair
            g_list[-1].add_edge(adding_ovlp[0], adding_ovlp[1], match=adding_ovlp[2])
            ovlp.pop(adding_ovlp)  # remove the added overlap pair
            ovlp.pop(rc_overlap(adding_ovlp))  # remove its rc too
            last_added_nodes = set([adding_ovlp[0], adding_ovlp[1]])
            if log: print('Added nodes this round: {0}'.format(last_added_nodes))
        if log: print('{0} graphs so far'.format(len(g_list)))

    # Assign subgraphs to DAGs and DCGs
    g_dict = {'dag':[], 'dcg':[]}
    for item in g_list:
        if nx.is_directed_acyclic_graph(item):
            g_dict['dag'].append(item)
        else:
            g_dict['dcg'].append(item)

    return g_dict


# trim the overlap_layout graph: Given a DiGraph, return a list of graph trimmed
# If the resulting graph is cyclic, turn it into acyclic
def trim_ovlp_graph(og):
    g = reduce_graph(og)  # og is a DiGraph() object, g is a DiGraph() object
    g = remove_branch_edges(g)
    g = subgraph(g)  # g became a list of DiGraph() objects
    g_dict = {'dag':[], 'dcg':[]}  # a dictionary with DAG and DCG for acyclic and cyclic graphs
    for item in g:
        if len(item.nodes) > 1:
            if nx.is_directed_acyclic_graph(item):  # if the subgraph is a DAG
                g_dict['dag'].append(item)
            else:  # Or otherwise it is a cycle or more
                item = subgraph(remove_branch_edges(item))  # trim the cycle first, item --> list
                for sg in item:
                    if nx.is_directed_acyclic_graph(sg):  # If trim change it to DAG
                        g_dict['dag'].append(reduce_graph(sg))
                    else:  # Still a cycle, return one cycle with longest path
                        cycles = list(nx.simple_cycles(sg))
                        cycles.sort(key=lambda z: len(z), reverse=True)
                        c1 = cycles[0]  # the longest cycle
                        x = nx.DiGraph()  # A new DAG for the broken cycle
                        for i in range(0, len(c1) - 1, 1):
                            x.add_edge(c1[i], c1[i + 1])
                        if len(x.nodes) > 1:
                            g_dict['dcg'].append(x)
    for item in g_dict['dag']:
        item.add_nodes_from(og.nodes(data=True))
        item.add_edges_from((u, v, og.edges[u, v]) for u, v in item.edges)
    for item in g_dict['dcg']:
        item.add_nodes_from(og.nodes(data=True))
        item.add_edges_from((u, v, og.edges[u, v]) for u, v in item.edges)
    return g_dict


# Return a list of node names to be removed from origin sequences
def rmlist_graph(g):
    rmlist = []
    for item in g:
        for node in item.nodes:
            rmlist.append(node[0])
    rmlist = set(rmlist)
    return rmlist


# Convert a list of graphs into a GFA format output
def graph_to_gfa(g_list, seqs):
    segment = []
    linker = []
    seq_len = {i[0]:len(i[1]) for i in seqIO.sequence(seqs)}
    for item in g_list:
        for node in item.nodes:
            segment.append(node[0])
        for edge in item.edges:
            output_line = ['L', edge[0][0], edge[0][1], edge[1][0], edge[1][1], str(item[edge[0]][edge[1]]['match'])]
            linker.append(output_line)
    segment = list(set(segment))
    segment = [['S', i, '*', 'LN:i:' + str(seq_len[i])] for i in segment]
    segment.sort()
    linker.sort()
    return segment + linker
