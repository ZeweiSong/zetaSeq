#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:36:47 2018

Bayesian multi-alignment basher

Coders who love to comment their code are unlikely to have bad luck.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
import networkx as nx

# Read in taxonomy file into a dictionary
def load_taxonomy(taxonomy_file):
    taxonomy = {}
    with open(taxonomy_file, 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            tid = line[0]
            taxon = tuple(line[1].split(';') + [tid])
            taxonomy[tid] = taxon
    return taxonomy

# Read in a blast6 format alignment
# Return a dictionary of targets
def load_blast6(alignment_file, taxonomy):
    target = {}
    with open(alignment_file, 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            qid = line[0]
            tid = taxonomy[line[1]]
            pident = float(line[2])
            try:
                target[tid].update({qid:pident})
            except KeyError:
                target[tid] = {qid:pident}
    return target

# Specify a taxonomic level to group all targets.
# If genus, all speceis will be grouped together.
def remove_redundant_target(tdict, threshold='species'):
    levels = {'kingdom':0, 'k':0, 'phylum':1, 'p':1, \
             'class':2, 'c':2, 'order':3, 'o':3, \
             'family':4, 'f':4, 'genus':5, 'g':5, \
             'species':6, 's':6, 'isolate':7, 'i':7}
    tout = {}
    try:
        trim_level = levels[threshold]
    except KeyError:
        print('{0} is not a eligible taxonomic level.'.format(threshold))
        print('Please choose from: {0}.',format(levels.keys()))
        return None
    for key, value in tdict.items():
        tid = tuple(key[:trim_level+1]) # get the new tid
        tout[tid] = tout.get(tid, {})
        for qid in value.keys(): # Update the query with the larger Pident value
            tout[tid][qid] = max(tout[tid].get(qid, 0), value[qid])
    return tout

# def target-query dictionary
def qtdict(alignment_file, taxonomy_file, threshold='species'):
    # read in taxonomy as a dictionary
    print('Read in taxonomy data from: {0}'.format(taxonomy_file))
    taxonomy = load_taxonomy(taxonomy_file)
    print('Found {0} record.'.format(len(taxonomy)))
    
    # Read in alignment, associate with taxonomy
    print('Read in alignment data from: {0}'.format(alignment_file))
    tdict = load_blast6(alignment_file, taxonomy)
    print('Found {0} targets in the alignment.'.format(len(tdict)))
    
    # Remove redundant target at the taxonomic level
    print('Compact the target to the level of {0}'.format(threshold))
    tout = remove_redundant_target(tdict, threshold=threshold)
    print('Targets were reduced from {0} to {1}'.format(len(tdict), len(tout)))
    
    return tout


# Build the query-target graph
def build_qtgraph(tdict):
    G = nx.Graph()
    tnames = {}
    qnames = {}
    for key, value in tdict.items():
        tnames[key] = tnames.get(key, 0) + len(value) # value is the number of linked queries
        for qid in value.keys():
            G.add_edge(key, qid, weight=value[qid])
            qnames[qid] = qnames.get(qid, 0) + 1 # value is the number of linked targets
    G.graph['tnames'] = tnames
    G.graph['qnames'] = qnames
    hit_dist = {}
    for key, value in qnames.items():
        hit_dist[value] = hit_dist.get(value, 0) + 1
    hits = sorted([(i,j) for i,j in hit_dist.items()], key=lambda x:x[0])
    print('Found total hits {0}'.format(sum([i[0]*i[1] for i in hits])))
    print('Query distribution (total = {0}):'.format(len(qnames)))
    print('\tHit\tCount')
    for item in hits:
        print('\t{0}\t{1}'.format(item[0], item[1]))
    return G


# Assign uniqueness to all target
# Use update_target_uniqueness() and remove_zero_uniqueness
def assign_uniqueness(qtg):
    qtg = update_target_uniqueness(qtg)
    print('Assigned uniqueness to {0} targets.'.format(len(qtg.graph['tnames'])))
    qtg = remove_zero_uniqueness(qtg)
    print('Kept {0} targets.'.format(len(qtg.graph['tnames'])))
    qtg = update_target_uniqueness(qtg)
    print('Update uniqueness to new set of targets.')
    return qtg

# Update the uniqueness for all targets
def update_target_uniqueness(qtg):
    for target in qtg.graph['tnames']:
        unique_count = 0
        for query in qtg.neighbors(target):
            if qtg.degree(query) == 1:
                unique_count += 1
            else:
                pass
        qtg.nodes[target]['uniqueness'] = unique_count
    return qtg

# Remove target with uniqueness == 0 (<1)
# After removing, some queris may have degree == 0
def remove_zero_uniqueness(qtg):
    remove_tnames = []
    for target in qtg.graph['tnames']:
        if qtg.nodes[target]['uniqueness'] < 1:
            remove_tnames.append(target)
            qtg.remove_node(target)
    for item in remove_tnames:
        del qtg.graph['tnames'][item]
    return qtg


# Update p(Ti|Qi) for all queries, after assigned uniqueness to all targets.
# p(Ti|Qi) = P(Qi|Ti)P(Ti)/P(Qi)
# assign new attribute 'probability' to queris
def update_query_probability(qtg):
    for query in qtg.graph['qnames']:
        if qtg.degree(query) > 1:
            pdict = {}
            for target in qtg.neighbors(query):
                pdict[target] = {'similarity':qtg.edges[query, target]['weight'],\
                                 'uniqueness':qtg.nodes[target]['uniqueness']}
            tsum = sum([i['uniqueness'] for i in pdict.values()])            
            for key, value in pdict.items():
                pdict[key]['product'] = value['similarity'] * (value['uniqueness']/tsum)
            psum = sum([i['product'] for i in pdict.values()])
            qtg.nodes[query]['probability'] = {i:j['product']/psum for i, j in pdict.items()} # add attribute to query
        elif qtg.degree(query) == 1:
            for target in qtg.neighbors(query):
                qtg.nodes[query]['probability'] = {target: 1}
    return qtg

# 
def qtgraph(tdict):
    pass

# Count the queirs into target profile
def count_profile(qtg):
    profile = {}
    for query in qtg.graph['qnames'].keys():
        try:
            for target, value in qtg.nodes[query]['probability'].items():
                profile[target] = profile.get(target, 0) + value
        except KeyError:
            pass
    return profile


#%% Test
align_file = 'alignment.txt'
taxon_file = 'taxonomy.txt'

target_dict = load_blast6(align_file, taxon_file)
#print(target_dict)
target_dict_nr = remove_redundant_target(target_dict, threshold='species')
#print(target_dict_nr)

qtGraph = build_qtgraph(target_dict_nr)
#print(qtGraph.nodes())
#%%
qtGraph = assign_uniqueness(qtGraph)
# for item in qtGraph.nodes():
#     print(item, qtGraph.nodes[item])
#%%
qtGraph = update_query_probability(qtGraph)
for item in qtGraph.graph['tnames']:
    print(item, qtGraph.nodes[item])
for item in qtGraph.graph['qnames']:
    print(item, qtGraph.nodes[item])
#%%
sample_profile = count_profile(qtGraph)
for key, value in sample_profile.items():
    print(key, value)