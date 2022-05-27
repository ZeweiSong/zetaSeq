# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Created on Thu May 19 17:33:42 2022

Dereplicate a group of contigs

dereco prep -i /path/to/contigs.fa.gz -o /path/to/out/folder/ -c 1000000
dereco cacl -i /path/to/out/folder -t 4

He who loves to comment his code is unlikely to have bad luck.
@author: Song, Zewei
@contact: songzewei@genomics.cn

"""

#%%
import random

# Generate a random sequence with length = l
def random_sequence(l):
    n = ('A', 'T', 'C', 'G')
    s = ''
    for i in range(l):
        s += random.sample(n, 1)[0]
    return s

def random_contig():
    content = []
    content.append(('t1', random_sequence(2000000))) # Add a 2M contig
    content += [('t2', content[0][1][3000:8000])] # t2 is a 5k subset of t1
    content += [('t3', content[1][1][1000:3500])] # t3 is a 2.5k subset of t2 aka t3 -> t2 -> t1
    content += [('t4', content[0][1][200000:202500])] # t4 = 2.5k
    content += [('t5', content[0][1][150000:208000])] # t4 -> t5 -> t1
    content.append(('t6', random_sequence(30000))) # add a 30k contig
    content += [('t7', content[5][1][1000: 5000])] # t7 -> t6
    content += [('t8', content[5][1][1500: 6500])] # t8 -> t6
    content.append(('t9',random_sequence(2100)))
    content.append(('t10', random_sequence(2100)))
    return content

#%%
for i in range(3):
    with open('tst_ctg_{0}.fa'.format(i+1), 'w') as f:
        ctg = random_contig()
        print('Generated tst_ctg_{0}.fa:'.format(i+1))
        for item in ctg:
            print('{0}: {1}'.format(item[0], len(item[1])))
            f.write('>{0}\n{1}\n'.format(item[0], item[1]))