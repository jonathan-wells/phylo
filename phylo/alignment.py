#!/usr/bin/env python3

import re
import numpy as np


class PairwiseAlignment():
    """Methods for pairwise alignments of sequences.

    Accepts protein nucleic acid sequence. Currently implemented:
        Needleman-Wunsch
        Smith-Waterman
        Needleman-Wunsch with affine gap penalties
        Smith-Waterman with affine gap penalties
    """

    def __init__(self, seq1, seq2, atype, scoremat):
        self.seq1 = seq1
        self.seq2 = seq2
        self.m = len(seq1)
        self.n = len(seq2)
        self.atype = atype

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

    def _clear_matrices(self):
        self.f_mat = np.zeros((self.m + 1, self.n + 1))
        self.ix_mat = np.zeros((self.m + 1, self.n + 1))
        self.iy_mat = np.zeros((self.m + 1, self.n + 1))

        self.f_ptr = np.zeros((self.m + 1, self.n + 1))
        self.ix_ptr = np.zeros((self.m + 1, self.n + 1))
        self.iy_ptr = np.zeros((self.m + 1, self.n + 1))

        self.ptr_mats = {'f': self.f_ptr, 'ix': self.ix_ptr, 'iy': self.iy_ptr}

    def _get_score(self, res1, res2):
        """Return score from specified scoring matrix"""
        return self.scoremat[self.resind[res1], self.resind[res2]]

    def _seq1_match(self, aln1, aln2, i, j):
        aln1 += self.seq1[i-1]
        aln2 += self.seq2[j-1]
        i -= 1
        j -= 1
        return aln1, aln2, i, j

    def _seq1_ins(self, aln1, aln2, i, j):
        aln1 += self.seq1[i-1]
        aln2 += '-'
        i -= 1
        return aln1, aln2, i, j

    def _seq1_del(self, aln1, aln2, i, j):
        aln1 += '-'
        aln2 += self.seq2[j-1]
        j -= 1
        return aln1, aln2, i, j

    def _set_affine_ptr(self, i, j,
                        fmatch, ixmatch, iymatch,
                        ixopen, ixextend,
                        iyopen, iyextend):
        if self.f_mat[i, j] == fmatch:
            self.f_ptr[i, j] = 2  # Go diagonal
        elif self.f_mat[i, j] == ixmatch:
            self.f_ptr[i, j] = 1  # Arrived from ix
        elif self.f_mat[i, j] == iymatch:
            self.f_ptr[i, j] = 3  # Arrived from iy
        if self.ix_mat[i, j] == ixopen:
            self.ix_ptr[i, j] = 4  # Arrived from f
        elif self.ix_mat[i, j] == ixextend:
            self.ix_ptr[i, j] = 5
        if self.iy_mat[i, j] == iyopen:
            self.iy_ptr[i, j] = 4  # Arrived from f
        elif self.iy_mat[i, j] == iyextend:
            self.iy_ptr[i, j] = 5

    def _traceback(self, aln1, aln2, i, j):
        """Decides necessary pointer in traceback matrix."""
        if self.f_ptr[i, j] == 2:
            step = self._seq1_match(aln1, aln2, i, j)
        elif self.f_ptr[i, j] == 1:
            step = self._seq1_ins(aln1, aln2, i, j)
        elif self.f_ptr[i, j] == 3:
            step = self._seq1_del(aln1, aln2, i, j)
        return step

    def _affine_traceback(self, currmat, aln1, aln2, i, j):
        """Decides necessary pointer in traceback matrices."""
        ptr_mat = self.ptr_mats[currmat]
        step = aln1, aln2, i, j
        if ptr_mat[i, j] == 2:
            step = self._seq1_match(*step)
        elif ptr_mat[i, j] == 1: # Go to ixtmat
            step = self._seq1_match(*step)
            currmat = 'ix'
        elif ptr_mat[i, j] == 3:
            step = self._seq1_match(*step)
            currmat = 'iy'
        elif ptr_mat[i, j] == 5 and currmat == 'ix':
            step = self._seq1_ins(*step)
        elif ptr_mat[i, j] == 5 and currmat == 'iy':
            step = self._seq1_del(*step)
        elif ptr_mat[i, j] == 4 and currmat == 'ix':
            step = self._seq1_ins(*step)
            currmat = 'f'
        elif ptr_mat[i, j] == 4 and currmat == 'iy':
            step = self._seq1_del(*step)
            currmat = 'f'
        aln1, aln2, i, j = step
        return currmat, aln1, aln2, i, j

    ############################################################################
    ## Algorithm implementations
    ############################################################################

    def nw_align(self, gap_open=8):
        """Needleman-Wunsch global alignment"""
        # Initialise matrices
        self._clear_matrices()
        for i in range(1, self.m + 1):
            self.f_mat[i, 0] = self.f_mat[i-1, 0] - gap_open
            self.f_ptr[i, 0] = 1
        for j in range(1, self.n+1):
            self.f_mat[0, j] = self.f_mat[0, j-1] - gap_open
            self.f_ptr[0, j] = 3
        # Populate scoring and traceback matrices
        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                a, b = self.seq1[i-1], self.seq2[j-1]
                match = self.f_mat[i-1, j-1] + self._get_score(a, b)
                seq1_ins = self.f_mat[i-1, j] - gap_open
                seq1_del = self.f_mat[i, j-1] - gap_open

                self.f_mat[i, j] = max(match, seq1_ins, seq1_del)
                if self.f_mat[i, j] == match:
                    self.f_ptr[i, j] = 2
                elif self.f_mat[i, j] == seq1_ins:
                    self.f_ptr[i, j] = 1
                elif self.f_mat[i, j] == seq1_del:
                    self.f_ptr[i, j] = 3

        # Traceback
        aln1, aln2 = '', ''
        i, j = self.m, self.n
        tracestep = aln1, aln2, i, j
        while i > 0 or j > 0:
            tracestep = self._traceback(*tracestep)
            i, j = tracestep[2:]
        aln1, aln2 = tracestep[:2]

        return aln1[::-1], aln2[::-1]

    def sw_align(self, gap_open=8):
        """Smith-Waterman local alignment"""
        # Initialise matrices
        self._clear_matrices()
        # Populate scoring and traceback matrices
        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                a, b = self.seq1[i-1], self.seq2[j-1]
                match = self.f_mat[i-1, j-1] + self._get_score(a, b)
                seq1_ins = self.f_mat[i-1, j] - gap_open
                seq1_del = self.f_mat[i, j-1] - gap_open

                self.f_mat[i, j] = max(match, seq1_ins, seq1_del, 0)
                if self.f_mat[i, j] == match:
                    self.f_ptr[i, j] = 2
                elif self.f_mat[i, j] == seq1_ins:
                    self.f_ptr[i, j] = 1
                elif self.f_mat[i, j] == seq1_del:
                    self.f_ptr[i, j] = 3

        # Traceback
        aln1, aln2 = '', ''
        i, j = np.unravel_index(np.argmax(self.f_mat, axis=None), self.f_mat.shape)
        tracestep = aln1, aln2, i, j
        while self.f_ptr[i, j] != 0:
            tracestep = self._traceback(*tracestep)
            i, j = tracestep[2:]
        aln1, aln2 = tracestep[:2]

        return aln1[::-1], aln2[::-1]

    def nwa_align(self, gap_open=8, gap_extend=1):
        """Needleman-Wunsch global alignment with affine gap penalty."""
        # Initialise matrices
        self._clear_matrices()
        for i in range(1, self.m + 1):
            self.f_mat[i, 0] = -gap_open - (i-1)*gap_extend
            self.ix_mat[i, 0] = -gap_open - (i-1)*gap_extend
            self.iy_mat[i, 0] = -gap_open - (i-1)*gap_extend
            self.f_ptr[i, 0] = 1
        for j in range(1, self.n + 1):
            self.f_mat[0, j] = -gap_open - (j-1)*gap_extend
            self.ix_mat[0, j] = -gap_open - (j-1)*gap_extend
            self.iy_mat[0, j] = -gap_open - (j-1)*gap_extend
            self.f_ptr[0, j] = 3

        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                a, b = self.seq1[i-1], self.seq2[j-1]

                ixopen = self.f_mat[i-1, j] - gap_open
                ixextend = self.ix_mat[i-1, j] - gap_extend
                self.ix_mat[i, j] = max(ixopen, ixextend)

                iyopen = self.f_mat[i, j-1] - gap_open
                iyextend = self.iy_mat[i, j-1] - gap_extend
                self.iy_mat[i, j] = max(iyopen, iyextend)

                fmatch = self.f_mat[i-1, j-1] + self._get_score(a, b)
                ixmatch = self.ix_mat[i-1, j-1] + self._get_score(a, b)
                iymatch = self.iy_mat[i-1, j-1] + self._get_score(a, b)
                self.f_mat[i, j] = max(fmatch, ixmatch, iymatch)

                self._set_affine_ptr(i, j,
                                     fmatch, ixmatch, iymatch,
                                     ixopen, ixextend,
                                     iyopen, iyextend)

        aln1, aln2 = '', ''
        i, j = self.m, self.n
        currmat = 'f'
        tracestep = currmat, aln1, aln2, i, j
        while i > 0 or j > 0:
            tracestep = self._affine_traceback(*tracestep)
            i, j = tracestep[-2:]
        aln1, aln2 = tracestep[1:3]

        return aln1[::-1], aln2[::-1]

    def swa_align(self, gap_open=8, gap_extend=1):
        self._clear_matrices()
        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                a, b = self.seq1[i-1], self.seq2[j-1]

                ixopen = self.f_mat[i-1, j] - gap_open
                ixextend = self.ix_mat[i-1, j] - gap_extend
                self.ix_mat[i, j] = max(ixopen, ixextend)

                iyopen = self.f_mat[i, j-1] - gap_open
                iyextend = self.iy_mat[i, j-1] - gap_extend
                self.iy_mat[i, j] = max(iyopen, iyextend)

                fmatch = self.f_mat[i-1, j-1] + self._get_score(a, b)
                ixmatch = self.ix_mat[i-1, j-1] + self._get_score(a, b)
                iymatch = self.iy_mat[i-1, j-1] + self._get_score(a, b)
                self.f_mat[i, j] = max(fmatch, ixmatch, iymatch, 0)

                self._set_affine_ptr(i, j,
                                     fmatch, ixmatch, iymatch,
                                     ixopen, ixextend,
                                     iyopen, iyextend)
        # Traceback
        currmat = 'f'
        aln1, aln2 = '', ''
        i, j = np.unravel_index(np.argmax(self.f_mat, axis=None), self.f_mat.shape)
        tracestep = currmat, aln1, aln2, i, j
        while self.f_ptr[i, j] != 0:
            tracestep = self._affine_traceback(*tracestep)
            aln1, aln2, i, j = tracestep[1:]
        return aln1[::-1], aln2[::-1]


def main():
    seq1 = 'ALFGLKSGRNGRITCMASYKVKLITPDGPIEFLFGLKSGRNGRITCMASYKVKLITPDGPECP'
    seq2 = 'ALFGLKLKRGDLAVAMASYKVDGTQEFECPLFGLKSGRNGRITCTCMASYKVKLITPDMASYKVKLITPDGP'
    test = PairwiseAlignment(seq1, seq2, 'prot', 'BLOSUM62')
    test.nw_align()
    test.nwa_align()
    test.sw_align()
    test.swa_align()

if __name__ == '__main__':
    main()
