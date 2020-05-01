#!/usr/bin/env python3

import numpy as np
from Bio import SeqIO

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
    trim5_seqs, trim5_msa = read_msa('../data/trim5.fna')
    print(trim5_seqs)
    
