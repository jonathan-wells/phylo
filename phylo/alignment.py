#!/usr/bin/env python3

import numpy as np
from Bio import SeqIO
import re


class PairwiseAlignment(object):
    """Methods for pairwise alignments of sequences.

    Accepts protein nucleic acid sequence. Currently implemented:
        Needleman-Wunsch
        Smith-Waterman
    """

    def __init__(self, seq1, seq2, atype, scoremat):
        self.seq1 = seq1
        self.seq2 = seq2
        self.m = len(seq1)
        self.n = len(seq2)
        
        self.scoremat = []
        self.resind = {}
        with open(f'../data/{scoremat}.txt') as infile:
            i = 0
            for line in infile:
                if re.match(r'^[A-z\*]', line):
                    line = line.split()
                    res, scores = line[0], [int(i) for i in line[1:]]
                    self.resind[res] = i
                    self.scoremat.append(scores)
                    i += 1
        self.scoremat = np.asarray(self.scoremat)

    ############################################################################ 
    ## Algorithm subroutines
    ############################################################################ 

    def get_score(self, res1, res2):
        """Return score from specified scoring matrix"""
        return self.scoremat[self.resind[res1], self.resind[res2]]

    def traceback(self, tmat, seq1_aln, seq2_aln, i, j):
        """Decides necessary pointer in traceback matrix."""
        if tmat[i, j] == 2:
            seq1_aln += self.seq1[i-1]
            seq2_aln += self.seq2[j-1]
            i -= 1
            j -= 1
        elif tmat[i, j] == 1:
            seq1_aln += self.seq1[i-1]
            seq2_aln += '-'
            i -= 1
        elif tmat[i, j] == 3:
            seq1_aln += '-'
            seq2_aln += self.seq2[j-1]
            j -= 1
        return seq1_aln, seq2_aln, i, j
    
    ############################################################################ 
    ## Algorithm implementations
    ############################################################################ 

    def nw_align(self, gap_open=8):
        """Needleman-Wunsch global alignment"""
        # Initialise matrices
        tmat = np.zeros((self.m + 1, self.n + 1))
        fmat = np.zeros((self.m + 1, self.n + 1))
        for i in range(self.m):
            fmat[i+1, 0] = fmat[i, 0] - gap_open
            tmat[i+1, 0] = 1
        for j in range(self.n):
            fmat[0, j+1] = fmat[0, j] - gap_open
            tmat[0, j+1] = 3
        
        # Populate scoring and traceback matrices
        for i in range(self.m):
            for j in range(self.n):
                a, b = self.seq1[i], self.seq2[j]
                match_mismatch = fmat[i, j] + self.get_score(a, b)
                seq1_ins = fmat[i, j+1] - gap_open
                seq1_del = fmat[i+1, j] - gap_open
                
                fmat[i+1, j+1] = max(match_mismatch, seq1_ins, seq1_del)
                if fmat[i+1, j+1] == match_mismatch:
                    tmat[i+1, j+1] = 2
                elif fmat[i+1, j+1] == seq1_ins:
                    tmat[i+1, j+1] = 1
                elif fmat[i+1, j+1] == seq1_del:
                    tmat[i+1, j+1] = 3

        # Traceback
        seq1_aln, seq2_aln = '', ''
        i, j = self.m, self.n
        tracestep = seq1_aln, seq2_aln, i, j
        while i > 0 or j > 0:
            tracestep = self.traceback(tmat, *tracestep)
            i, j = tracestep[2:]
        seq1_aln, seq2_aln = tracestep[:2]
        
        return seq1_aln[::-1], seq2_aln[::-1]

    def sw_align(self, gap_open=8):
        """Smith-Waterman local alignment"""
        # Initialise matrices
        tmat = np.zeros((self.m + 1, self.n + 1))
        fmat = np.zeros((self.m + 1, self.n + 1))
        
        # Populate scoring and traceback matrices
        for i in range(self.m):
            for j in range(self.n):
                a, b = self.seq1[i], self.seq2[j]
                match_mismatch = fmat[i, j] + self.get_score(a, b)
                seq1_ins = fmat[i, j+1] - gap_open
                seq1_del = fmat[i+1, j] - gap_open
                
                fmat[i+1, j+1] = max(match_mismatch, seq1_ins, seq1_del, 0)
                if fmat[i+1, j+1] == match_mismatch:
                    tmat[i+1, j+1] = 2
                elif fmat[i+1, j+1] == seq1_ins:
                    tmat[i+1, j+1] = 1
                elif fmat[i+1, j+1] == seq1_del:
                    tmat[i+1, j+1] = 3

        # Traceback
        seq1_aln, seq2_aln = '', ''
        i, j = np.unravel_index(np.argmax(fmat, axis=None), fmat.shape)
        tracestep = seq1_aln, seq2_aln, i, j
        while tmat[i, j] != 0:
            tracestep = self.traceback(tmat, *tracestep)
            i, j = tracestep[2:]
        seq1_aln, seq2_aln = tracestep[:2]

        return seq1_aln[::-1], seq2_aln[::-1]

    def swa_align(self, gap_open=8, gap_extend=1):
        pass

def read_msa(filename):
    seqids = []
    array = []
    for record in SeqIO.parse(filename, 'fasta'):
        seqids.append(record.name)
        array.append(str(record.seq).split())
    return seqids, np.asarray(array)

def calc_distance(seq1, seq2):
    pass

if __name__ == '__main__':
    test = PairwiseAlignment('HEAGAWGHEE',
                             'PAWHEAE',
                             'prot',
                             'BLOSUM50')
    s1, s2 = test.nw_align()
    print(s1)
    print(s2)
    s1, s2 = test.sw_align()
    print(s1)
    print(s2)
