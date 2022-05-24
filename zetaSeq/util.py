#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:36:47 2018

Coders who love to comment their code are unlikely to have bad luck.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%
from __future__ import print_function
from __future__ import division

def test_me(args):
    from zetaSeq import io as seqIO
    seqIO.test_me()
    print('Inputfile assigned as {0}'.format(args.inputfile))

def concat_biom():
    pass


#%% count the number of reads aligned to the reference.
# Reads may align to more than one record.
def countAlnReads(args):
    inputFile = args.input
    
    reads = {}
    with open(inputFile, 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            currentRead = line[0]
            reads[currentRead] = reads.get(currentRead, 0) + 1    
    
    print('There are {0} reads in the alignment {1}.'.format(len(reads), inputFile))
    return None

#%% keep paired alignments
def keepPairAln(args):
    outputAln = args.output
    inputAln = []
    for f in args.input:
        inputAln.append(f)
    if len(inputAln) != 2:
        print('The number of files is not equal to 2.')
        return None
    else:
        pass
    # Build the query dictionary
    # Two alignments are combined into one
    
    queryDict = {}
    print('Parsing Read1 alignment ...')
    with open(inputAln[0], 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            line[0] = line[0].split('/')[0] # Get the Query label (stripped Read info)
            try:
                queryDict[line[0]]['r1'].append(line[1])
            except KeyError:
                queryDict[line[0]] = {'r1':[], 'r2':[]}
                queryDict[line[0]]['r1'].append(line[1])
    print('Found {0} queries in Read1 alignment.'.format(len(queryDict)))
    print('Parsing Read2 alignment ...')

    r2notinr1 = 0
    with open(inputAln[1], 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            line[0] = line[0].split('/')[0]
            try:
                queryDict[line[0]]['r2'].append(line[1])
            except KeyError:
                queryDict[line[0]] = {'r1':[], 'r2':[]}
                queryDict[line[0]]['r2'].append(line[1])
                r2notinr1 += 1
    print('Total queries add to {0}'.format(len(queryDict)))
    
    aln = {}
    for query, alignments in queryDict.items():
        try:
            overlap = tuple(set(alignments['r1']).intersection(set(alignments['r2'])))
            if len(overlap) > 0:
                aln[query] = overlap
        except KeyError:
            pass
    
    print('Found {0} queries in both alignment.'.format(len(aln)))
    print('\t representing {0} alignments.'.format(sum([len(i) for i in aln.values()])))
    print('Found {0} queries ONLY in Read2 alignment.'.format(r2notinr1))
    
    with open(outputAln, 'w') as f:
        for query, alignments in aln.items():
            for item in alignments:
                f.write('{0}\t{1}\n'.format(query, item))
    print('Paired alignment wrote to {0}.'.format(outputAln))
    
    return None


#%% concat multiple biome files into one biom file, using JSON format
def concatBioms(args):
    from biom.table import Table
    import json
    
    biomFile = args.output
    biomList = []
    for f in args.input:
        biomList.append(f)
    if len(biomList) <= 1:
        print('Found only one biom profile, will still give you a (brand) new one.')
        biomProfile = Table.from_json(json.load(biomList[0]))
        with open(biomFile, 'w') as f:
            biomProfile.to_json('Generated_by_almighty_metaSeq', f)
    else:
        print('Found {0} biom profiles under.'.format(len(biomList)))
        biomProfile = Table.from_json(json.load(biomList[0]))
    
        for f in biomList[1:]:
            biomProfile = biomProfile.concat([Table.from_json(json.load(f))])
    
        with open(biomFile, 'w') as f:
            biomProfile = biomProfile.sort()
            biomProfile.to_json('Generated_by_almighty_metaSeq', f)
        print('Concatenated {0} profiles into {1}.'.format(len(biomList), biomFile))    
    return None

#%% Winner take all profile
def winnerTakeAll(args):
    from zetaSeq import amplicon
    alnString = []
    for item in args.input:
        alnString.append(item)
    # targetNumber = len(alnString) # The number of targets used (normally this is 2, e.g. ITS1 forward & ITS1 reverse)
    sampleName = args.sample_name
    tabFile = args.output
    
    greedy = args.greedy
    less_greedy = args.less_greedy
    if less_greedy:
        greedy = False
    
    ec = args.ec
    median = args.median
    ave = args.ave
    if median:
        weight = 'median'
    elif ave:
        weight = 'average'
    elif ec:
        weight = 'ec'
    
    print('Ah the bloody capitalism came to {0}!'.format(sampleName))
    mode = 'greedy'
    if not greedy:
        mode = 'less greedy'
    print('Released are the {0} capitalists.'.format(mode))
    if not greedy:
        print("\tI'll spit out some profit to maximize mine.\n")
    else:
        print('\tAll profit is mine!\n')
    alnNormalized = amplicon.initAlignment(alnString)
    wtaProfile = amplicon.winnerTakeAll(alnNormalized, progress=True, greedy=greedy, weight=weight)
    print('Winner take all profile found! {0} references survived.'.format(len(wtaProfile)))
    profile = list(wtaProfile.items())
    profile.sort(key=lambda x:x[1], reverse=True)
    
    print('Tab-delimited profile wrote to {0}.'.format(tabFile))
    with open(tabFile, 'w') as f:
        f.write('{0}\t{1}\n'.format('Reference', sampleName))
        for item in profile:
            if int(item[1]) >= 1:
                f.write('{0}\t{1}\n'.format(item[0], int(item[1])))
    
    return None


#%% Randomly sample a profile array into lower depth
def rarefy(args):
    from zetaSeq import rarefy as rs
    from zetaSeq import io as seqIO
    import time
    
    input_otu = args.input
    
    iter_num = int(args.iteration)
    thread = int(args.threads)
    meta_col = args.meta_column
    output_otu = args.output

    otu_table = seqIO.parser_otu_table(input_otu, meta_col=meta_col)
    input_sample = otu_table.sample_matrix
    start = time.time()

    print('Input OTU table: %s' % input_otu)
    if args.depth:
        depth = int(args.depth)
        print('Sampling depth: %i' % depth)
    else:
        depth = min([sum(i[1:]) for i in input_sample])
        print('Sampling depth set to the minimum abundance: %i' % depth)
    print('Iteration time for each OTU: %i' % iter_num)
    print('Threads number: %i' % thread)
    print('Reading in the OTU table ...')

    if args.keep_all:
        count = 0
        for line in input_sample:
            if sum(line[1:]) < depth:
                count += 1
        print('Found %i samples in the OTU table.' % len(input_sample))
        print('%i samples has total abundance less than the sampling depth, but will be kept in the output.' % count)
    else:
        temp = []
        count = 0
        for line in input_sample:
            if sum(line[1:]) >= depth:
                temp.append(line)
            else:
                count += 1
        input_sample = temp
        print('Found %i samples in the OTU table.' % len(input_sample))
        print('%i samples has total abundance less than the sampling depth, and will be excluded.' % (count))

    otu_id = otu_table.species_id
    otu_table_rarefied = [['OTU_ID'] + otu_id + otu_table.meta_id]

    for sample in input_sample:
        print('Rarefying %s ...' % sample[0])
        repeat_sample = rs.repeat_rarefaction_parallel(sample[1:], depth, iter_num, processor=thread)
        repeat_sample.sort(key=lambda x: sum(i > 0 for i in x))
        repeat_sample = [sample[0]] + repeat_sample[int(iter_num / 2)]  # Pick the rarefied sample with the average richness
        otu_table_rarefied.append(repeat_sample[:])

    otu_table_rarefied = [list(i) for i in zip(*otu_table_rarefied)]

    # Add meta data
    meta_data = otu_table.meta_dict()
    for line in otu_table_rarefied[1:]:
        for key in otu_table.meta_id:
            line.append(meta_data[key][line[0]])
    for key in otu_table.meta_id:
        otu_table_rarefied[0].append(key)

    with open(output_otu, 'w') as f:
        for line in otu_table_rarefied:
            line = [str(i) for i in line]
            f.write('%s\n' % '\t'.join(line))
    print('Rarefied OTU table saved in %s.' % output_otu)
    end = time.time()
    used_time = round(float(end - start), 2)
    time_per_sample = round(used_time / len(input_sample), 2)
    print('Total time used: %s seconds (%s seconds per sample)' % (str(used_time), str(time_per_sample)))    
    
    return None