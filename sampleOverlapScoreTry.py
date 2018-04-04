#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 14:26:16 2018
@author: samhitapn
@purpose: Overlap score trials
@modified: 29.March.2018
"""

#import biopython
from Bio import pairwise2 
from Bio import SeqIO as sio
import numpy
import itertools

# Read in the fastq files
"""
for s in sio.parse("ONT_Sample1_5Reads.fastq","fastq"):
    print(s.seq[1])
    print(s.letter_annotations["phred_quality"][1])
    print("---------------- \n")
    #print(s.format("fastq"))
"""
seq = sio.to_dict(sio.parse("data/ONT_Sample1_5Reads.fastq","fastq"))
"""
for seq1 in seq:
    for seq2 in seq:
        if seq1 != seq2:
"""
for seq1, seq2 in itertools.combinations(seq, 2):
    al = pairwise2.align.globalxx(seq[seq1].seq,seq[seq2].seq)
    print(al[1])
    #print(pairwise2.format_alignment(*al[0]))
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
    



#alignments = pairwise2.align.globalxx("ACCGT", "ACG")
#print(alignments[1])
#print(pairwise2.format_alignment(*alignments))


