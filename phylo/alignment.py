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
        for i in range(1, self.m + 1):
            fmat[i, 0] = fmat[i-1, 0] - gap_open
            tmat[i, 0] = 1
        for j in range(self.n):
            fmat[0, j] = fmat[0, j-1] - gap_open
            tmat[0, j] = 3
        
        # Populate scoring and traceback matrices
        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                a, b = self.seq1[i-1], self.seq2[j-1]
                match = fmat[i-1, j-1] + self.get_score(a, b)
                seq1_ins = fmat[i-1, j] - gap_open
                seq1_del = fmat[i, j-1] - gap_open
                
                fmat[i, j] = max(match, seq1_ins, seq1_del)
                if fmat[i, j] == match:
                    tmat[i, j] = 2
                elif fmat[i, j] == seq1_ins:
                    tmat[i, j] = 1
                elif fmat[i, j] == seq1_del:
                    tmat[i, j] = 3

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
        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                a, b = self.seq1[i-1], self.seq2[j-1]
                match = fmat[i-1, j-1] + self.get_score(a, b)
                seq1_ins = fmat[i-1, j] - gap_open
                seq1_del = fmat[i, j-1] - gap_open
                
                fmat[i, j] = max(match, seq1_ins, seq1_del, 0)
                if fmat[i, j] == match:
                    tmat[i, j] = 2
                elif fmat[i, j] == seq1_ins:
                    tmat[i, j] = 1
                elif fmat[i, j] == seq1_del:
                    tmat[i, j] = 3

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
        tmat = np.zeros((self.m + 1, self.n + 1))
        ixtmat = np.zeros((self.m + 1, self.n + 1))
        iytmat = np.zeros((self.m + 1, self.n + 1))
        
        fmat = np.zeros((self.m + 1, self.n + 1))
        ixmat = np.zeros((self.m + 1, self.n + 1))
        iymat = np.zeros((self.m + 1, self.n + 1))
        
        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                a, b = self.seq1[i-1], self.seq2[j-1]

                ixopen = fmat[i-1, j] - gap_open
                ixextend = ixmat[i-1, j] - gap_extend
                ixmat[i, j] = max(ixopen, ixextend)
                
                iyopen = fmat[i, j-1] - gap_open
                iyextend = iymat[i, j-1] - gap_extend
                iymat[i, j] = max(iyopen, iyextend)

                fmatch = fmat[i-1, j-1] + self.get_score(a, b)
                ixmatch = ixmat[i-1, j-1] + self.get_score(a, b)
                iymatch = iymat[i-1, j-1] + self.get_score(a, b)
                fmat[i, j] = max(fmatch, ixmatch, iymatch, 0)
                
                if fmat[i, j] == fmatch:
                    tmat[i, j] = 2  # Go diagonal
                elif fmat[i, j] == ixmatch:
                    tmat[i, j] = 1  # Arrived from ix
                elif fmat[i, j] == iymatch:
                    tmat[i, j] = 3  # Arrived from iy
                if ixmat[i, j] == ixopen:
                    ixtmat[i, j] = 4  # Arrived from f
                elif ixmat[i, j] == ixextend:
                    ixtmat[i, j] = 5
                if iymat[i, j] == iyopen:
                    iytmat[i, j] = 4  # Arrived from f
                elif iymat[i, j] == iyextend:
                    iytmat[i, j] = 5
        # Traceback

        seq1_aln, seq2_aln = '', ''
        i, j = np.unravel_index(np.argmax(fmat, axis=None), fmat.shape)
        currmat = tmat
        cname = 'tmat'
        while currmat[i, j] != 0:
            if currmat[i, j] == 2:
                seq1_aln += self.seq1[i-1]
                seq2_aln += self.seq2[j-1]
                i -= 1
                j -= 1
            elif currmat[i, j] == 1: # Go to ixtmat
                currmat = ixtmat
                cname = 'ixtmat'
                seq1_aln += self.seq1[i-1]
                seq2_aln += self.seq2[j-1]
                i -= 1
                j -= 1
            elif currmat[i, j] == 3:
                currmat = iytmat
                cname = 'iytmat'
                seq1_aln += self.seq1[i-1]
                seq2_aln += self.seq2[j-1]
                i -= 1
                j -= 1
            elif currmat[i, j] == 5 and cname == 'ixtmat':
                seq1_aln += self.seq1[i-1]
                seq2_aln += '-'
                i -= 1
            elif currmat[i, j] == 5 and cname == 'iytmat':
                seq1_aln += '-'
                seq2_aln += self.seq2[j-1]
                j -= 1
            elif currmat[i, j] == 4 and cname == 'ixtmat':
                currmat = tmat
                cname = 'tname'
                seq1_aln += self.seq1[i-1]
                seq2_aln += '-'
                i -= 1
            elif currmat[i, j] == 4 and cname == 'iytmat':
                currmat = tmat
                cname = 'tname'
                seq1_aln += '-'
                seq2_aln += self.seq2[j-1]
                j -= 1
            
        return seq1_aln[::-1], seq2_aln[::-1]

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
    a1 = 'ALFGLKSG--RNGRITCMASYKVKLITPDGPIEFECP'
    a2 = 'GFLGLKTSLKRGDLAVAMASYKV-----DGTQEFECP'
    print(a1)
    print(a2)
    print()
    s1 = 'ALFGLKSGRNGRITCMASYKVKLITPDGPIEFECP'
    s2 = 'GFLGLKTSLKRGDLAVAMASYKVDGTQEFECP'
    test = PairwiseAlignment(s1, s2, 'prot', 'BLOSUM62')
    a1, a2 = test.swa_align()
    print(a1)
    print(a2)
    print()
    a1, a2 = test.sw_align()
    print(a1)
    print(a2)
    a1, a2 = test.nw_align()
    print(a1)
    print(a2)
