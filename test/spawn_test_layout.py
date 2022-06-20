#!/usr/bin/env python
# -*- coding: utf-8 -*-

from zetaSeq import io as seqIO
import random

# Generate a random sequence with length = l
def random_sequence(l):
    n = ('A', 'T', 'C', 'G')
    s = ''
    for i in range(l):
        s += random.sample(n, 1)[0]
    return s

# Reverse compliment a sequence
def rc(seq):
    d = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    seq_rc = ''.join([d[i] for i in seq])[::-1]
    return seq_rc

s = random_sequence(100000)
s0 = ('s0', s)
s1 = ('s1', rc(s[0:5000]))
s2 = ('s2', s[3000:8000])
s3 = ('s3', s[4000:12000])
s4 = ('s4', s[8000:18000])
s5 = ('s5', rc(s[14000:25000]))
s6 = ('s6', s[20000:30000])
s7 = ('s7', s[26000:40000])
s8 = ('s8', s[50000:55000])
s9 = ('s9', s[60000:70000])
s10 = ('s10', s[80000:90000])
ss = [s1, s2, s3, s4, s5, s6, s7, s8, s9, s10]
c = seqIO.write_seqs(ss, 'sim.fa', gz=False)
# c = seqIO.write_seqs([s0], 'ref.fa', gz=False)
