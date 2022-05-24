#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 16:59:32 2017

This script contains functions for handling taxonomy.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division


# Given a tax string, check Unclassified (and Candida) redundancy
# Return the fixed string
# The fixing principle is based on NCBI taxonomy dump
def fix_unclassified(taxa):
    if type(taxa) != 'list':
        taxa = taxa.split(';')
    assert len(taxa) == 7
    name = taxa[0][3:]
    taxa_nr = []
    for i, item in enumerate(taxa): # Starting from Phylum
        if i == 0:
            taxa_nr.append(item)
        else:
            if item[3:] == 'Unclassified' and taxa[i-1][3:15] != 'Unclassified':
                taxa_nr.append(item + taxa[i-1][3:])
                name = taxa[i-1][3:]
            elif item[3:] == 'Unclassified' and taxa[i-1][3:] == 'Unclassified':
                taxa_nr.append(item + name)
            elif item[3:] == 'Candida':
                taxa_nr.append(item + taxa[i-1][3:])
            else:
                taxa_nr.append(item)
    if taxa == taxa_nr:
        value = 0
    else:
        value = 1
    return (taxa_nr, value)


# Check the redundancy of a list of taxonomy strings, and return a report
# The taxa string has to be in seven levels (as a list or tuple)
def check_redundancy(tax_strings):
    length = {}
    for item in tax_strings:
        length[len(item)] = length.get(len(item), 0) + 1
    if len(length) == 1:
        taxa_len = tuple(length.values())[0]
    else:
        return length
    
    for i in range(taxa_len - 1):
        level = i + 1
        d = {}
        for item in tax_strings:
            d[item[level]][item[:level]] = d[item[level]].get(item[:level], 0) + 1
        redund = {}
        for key, value in d.items():
            if len(value) > 1:
                redund[key] = value
        if len(redund) > 0:
            return redund
            
    return True