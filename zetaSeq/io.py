#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 16:59:32 2017

This script contains functions that take a sequence file as input.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
from __future__ import print_function
from __future__ import division
import json

def test_me():
    print('This is a test for zetaSeq.')
    return None


# Function for checking the file type (gz, fasta, fastq)
# Return a tuple with two boolean value (T/F, T/F) --> (gz/not_gz, FASTQ/FASTA)
# FBI WARNING, this is a weak checker
    # It only check the file extension for gz file
    # It only check the first character for sequence type
def showMeTheType(filePath):
    import gzip
    seqType = {'>':False, '@':True}
    fileType = []
    if filePath.endswith('.gz'):
        fileType.append(True)
    else:
        fileType.append(False)
    if fileType[0]:
        with gzip.open(filePath, 'rt') as f:
            headSymbol = f.readline()[0]
    else:
        with open(filePath, 'r') as f:
            headSymbol = f.readline()[0]
    fileType.append(seqType[headSymbol])
    return tuple(fileType)


# an iterator oobject for reading a single sequence file. This is the most common Class.
# It will NOT check the format of the file, either can it deal with multiple line FASTA file.
# At my desktop computer, 1 M reads can be read in in about 2 seconds.
class sequence(object):
    def __init__(self, filePath):
        fileType = showMeTheType(filePath)
        if fileType[0]:
            import gzip
            self.file = gzip.open(filePath, 'rt')
        else:
            self.file = open(filePath, 'r')

        if fileType[1]:
            self.n = 4
        else:
            self.n = 2

    def __iter__(self):
        return self

    def __next__(self):
        record = []
        for i in range(self.n):
            line = self.file.readline().strip('\n')
            if line:
                record.append(line)
            else:
                raise StopIteration
        record[0] = record[0][1:]
        return record

# This is for FASTA with multilpe lines
        # Sequences are read in as a tuple of tuples
def sequence_multiline_fasta(filePath):
    fileType = showMeTheType(filePath)
    if fileType[0]:
        import gzip
        file = gzip.open(filePath, 'rt')
    else:
        file = open(filePath, 'r')    
    
    content = file.readlines() # Read in all lines
    content = [i.strip('\n') for i in content] # Remove EOF symbol
    fasta = []
    pos = [index for index, item in enumerate(content) if item[0] == ">"] # Get the position of all headers
    fasta = [(content[i][1:], ''.join(content[i+1:pos[index+1]])) for index, i in enumerate(pos[:-1])] # Get the sequence except for the last one
    fasta.append((content[pos[-1]][1:], ''.join(content[pos[-1]+1:]))) # Add the last one
    file.close()
    return fasta


# Same iterator but read in multiple record into memory at once
# By testing, read in file in trunk does not boost the reading speed.
class sequence_trunk(object):
    def __init__(self, filePath, trunk_size=100000):
        fileType = showMeTheType(filePath)
        if fileType[0]:
            import gzip
            self.file = gzip.open(filePath, 'rt')
        else:
            self.file = open(filePath, 'r')

        if fileType[1]:
            self.n = 4
        else:
            self.n = 2
        self.trunk_size = trunk_size

    def __iter__(self):
        return self

    def __next__(self):
        record_trunk = []
        for record in range(self.trunk_size):
            record = []
            for i in range(self.n):
                line = self.file.readline().strip('\n')
                if line:
                    record.append(line)
                else:
                    if len(record_trunk) > 0:
                        return record_trunk
                    else:
                        raise StopIteration
            record[0] = record[0][1:]
            record_trunk.append(record)
        return record_trunk


# Iterator for two files
# It only work for files with ABSOLUTELY corresponding record.
class sequence_twin(object):
    def __init__(self, file_r1, file_r2):
        fileType = showMeTheType(file_r1)
        fileType_2 = showMeTheType(file_r2)
        if fileType[0] != fileType_2[0] or fileType[1] != fileType_2[1]:
            print('Inconsistent file type, are you serious?')
            self.n = 1
        if fileType[0]:
            import gzip
            self.r1 = gzip.open(file_r1, 'rt')
            self.r2 = gzip.open(file_r2, 'rt')
        else:
            self.r1 = open(file_r1, 'r')
            self.r2 = open(file_r2, 'r')
        if fileType[1]: self.n = 4
        else: self.n = 2

    def __iter__(self):
        return self

    def __next__(self):
        record = [[],[]]
        for i in range(self.n):
            line_r1 = self.r1.readline().strip('\n')
            line_r2 = self.r2.readline().strip('\n')
            if line_r1:
                record[0].append(line_r1)
                record[1].append(line_r2)
            else:
                raise StopIteration
        record[0][0] = record[0][0][1:]
        record[1][0] = record[1][0][1:]
        return record[0], record[1]


class sequence_twin_trunk(object):
    def __init__(self, file_r1, file_r2, trunk_size=100000):
        fileType = showMeTheType(file_r1)
        fileType_2 = showMeTheType(file_r2)
        if fileType[0] != fileType_2[0] or fileType[1] != fileType_2[1]:
            print('Inconsistent file type, are you serious?')
            self.n = 1
        if fileType[0]:
            import gzip
            self.r1 = gzip.open(file_r1, 'rt')
            self.r2 = gzip.open(file_r2, 'rt')
        else:
            self.r1 = open(file_r1, 'r')
            self.r2 = open(file_r2, 'r')
        if fileType[1]: self.n = 4
        else: self.n = 2
        self.trunk_size = trunk_size

    def __iter__(self):
        return self

    def __next__(self):
        r1_trunk = []
        r2_trunk = []
        for record in range(self.trunk_size):
            r1 = []
            r2 = []
            for i in range(self.n):
                line_r1 = self.r1.readline().strip('\n')
                line_r2 = self.r2.readline().strip('\n')
                if line_r1:
                    r1.append(line_r1)
                    r2.append(line_r2)
                else:
                    if len(r1_trunk) > 0:
                        return r1_trunk, r2_trunk
                    else:
                        raise StopIteration
            r1[0] = r1[0][1:]
            r2[0] = r2[0][1:]
            r1_trunk.append(r1)
            r2_trunk.append(r2)
        return r1_trunk, r2_trunk


# This function is still under test, it reads in file bytes, should be a bit faster than read in by line
# Need to find a way to remove header symbol
# This is 40% faster than reading by line.
class sequence_bytes(object):
    def __init__(self, filePath, size = 1000000,fastx='a'):
        self.file = open(filePath, 'r')
        self.size = size
        self.fastx = fastx
        self.tail = ''

        if fastx == 'a':
            self.n = 2
        elif fastx == 'q':
            self.n = 4
        else:
            print('Please specify the right format, "a" for FASTA and "q" for FASTQ.')
    def __iter__(self):
        return self

    def __next__(self):
        content = self.file.read(self.size)

        if content:
            content = self.tail + content
            self.tail = ''
            content = content.split('\n')
            self.tail = content[-1]
            tail_n = (len(content) - 1) % 4 # Set the line step to 4, set to 2 for FASTA
            if tail_n > 0:
                self.tail = '\n'.join(content[:-1][-1*tail_n::]) + '\n' + self.tail
            else:
                self.tail = self.tail
            #print(self.tail)
            content = content[:(len(content) - tail_n - 1)]
            content = [content[x:x+self.n] for x in range(0, len(content), self.n)]
            return content
        else:
            if self.tail:
                content = self.tail.split('\n')
                content = [content[x:x+self.n] for x in range(0, len(content), self.n)]
                return content
            else:
                raise StopIteration


# Write the content to a fastx file
def write_seqs(seq_content, filePath, fastx='a', mode='w', gz=True):
    import gzip
    count = 0
    if fastx == 'a':
        header = '>'
    elif fastx == 'q':
        header = '@'
    else:
        header = '-_-b'
    if gz:
        if mode == 'w': mode = 'wt'
        else: mode = 'at'
        f = gzip.open(filePath, mode, newline='')
    else:
        f = open(filePath, mode, newline='')
    for record in seq_content:
        label = header + record[0]
        for line in [label] + list(record[1:]):
            f.write('%s\n' % line)
        count += 1
    f.close()
    return count


#%% Functions for alignment
class alignment(object):
    def __init__(self, alnFile):
        self.fileName = alnFile
        self.aln = open(alnFile, 'r')

    def __iter__(self):
        return self

    def __next__(self):
        line = self.aln.readline()
        if line:
            line = line.strip('\n').split('\t')
            return line
        else:
            raise StopIteration

def readAlignment(alnFile, sort=True):
    aln = alignment(alnFile)
    alnList = []
    for line in aln:
        alnList.append(line)
    if sort: # Sort using the first column, query by default.
        alnList.sort(key=lambda x:x[0])
    return alnList

def writeAlignment(alnList, outFile):
    i = 0
    with open(outFile, 'w') as f:
        for line in alnList:
            f.write('{0}\n'.format('\t'.join([str(i) for i in line])))
            i += 1
    return i


#%%
# Return reverse compliment of a sequence
# This part is got from Stakoverflow
#(https://stackoverflow.com/questions/19570800/reverse-complement-dna) by corinna
# Works for Python 3
def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhdN', 'TGCAtgcaYRKMyrkmBVDHbvdhN'))[::-1]


# Count the line number of a file
# Count the line breaker \n by reading a segment in bytes
# The code is borrowed from Stackoverflow
def blocks(inputFile, size=65536):
    while True:
        b = inputFile.read(size)
        if not b: break
        yield b

def count_line(inputFile):
    with open(inputFile, "r", encoding="utf-8", errors='ignore') as f:
        return sum(bl.count("\n") for bl in blocks(f))


#%% Functions for OTU table
class parser_otu_table(object):
    def __init__(self, file_path, meta_col='taxonomy'):
        with open(file_path, 'r') as f:
            temp = f.readlines()
        table = []
        for line in temp:
            table.append(line.strip('\n').split('\t'))
        size = len(table[0])
        for line in table:
            if len(line) < size:
                line += [''] * (size - len(line))

        # Get id for sample, species, amd meta data
        try:
            meta_position = table[0].index(meta_col)  # Begining position of meta data
        except ValueError:
            meta_position = len(table[0]) + 1
        self.sample_id = table[0][1:meta_position]  # Second column till the first meta column
        self.meta_id = table[0][meta_position:]  # Start from the first meta column to the end
        self.species_id = [i[0] for i in table[1:]]  # First column starting from the second row
        self.header = ['OTU_ID'] + self.sample_id + self.meta_id

        # Convert all abundance to intger
        for line in table[1:]:
            try:
                line[1:meta_position] = map(int, line[1:meta_position])
            except ValueError:
                print("There are non-number value in your OTU table.")
                import sys
                sys.exit()

        # Get sample, species, and meta data matrix with head names
        self.species_matrix = [i[:meta_position] for i in table[1:]]
        temp = [i[1:meta_position] for i in table]
        self.sample_matrix = [list(i) for i in zip(*temp)]
        temp = [i[meta_position:] for i in table]
        self.meta_matrix = [list(i) for i in zip(*temp)]


    # Generate a dictionary using sample name, OTU name
    def sample_dict(self):
        sample = {}
        for s in self.sample_id:
            sample[s] = {}
        for line in self.sample_matrix:
            abundance = line[1:]
            for i in range(len(self.species_id)):
                if int(abundance[i]) > 0:
                    sample[line[0]][self.species_id[i]] = int(abundance[i])
        return sample

    # Generate a dictionary using OTU name, sample name
    def species_dict(self):
        species = {}
        for s in self.species_id:
            species[s] = {}
        for line in self.species_matrix:
            for i in range(len(self.sample_id)):
                if int(line[1:][i]) > 0:
                    species[line[0]][self.sample_id[i]] = int(line[1:][i])
        return species


    # Generate a dictionary using mata name, OTU name
    def meta_dict(self):
        meta = {}
        for m in self.meta_id:
            meta[m] = {}
        for line in self.meta_matrix:
            for i in range(len(self.species_id)):
                meta[line[0]][self.species_id[i]] = line[1:][i]
        return meta
