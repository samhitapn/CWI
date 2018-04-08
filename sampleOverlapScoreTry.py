#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 14:26:16 2018
@author: samhitapn
@purpose: Overlap score trials
@modified: 04.April.2018
"""

#import biopython
from Bio import pairwise2 
from Bio import SeqIO as sio
import numpy
import itertools

# Read in the fastq files -> Usually the overalp pairs; here only the test 5 reads
sequences = sio.to_dict(sio.parse("data/ONT_Sample1_5Reads.fastq","fastq"))

# Get the alignment of all the overlap pairs
alignment = {}
pairCount = 1
for seq1, seq2 in itertools.combinations(sequences, 2):
    al = pairwise2.align.globalxx(sequences[seq1].seq,sequences[seq2].seq)
    alignment[pairCount] = sequences[seq1].id + "," + al[0][0] + "," + str(sequences[seq1].letter_annotations["phred_quality"]) + "," + sequences[seq2].id + "," + al[0][1] +  "," + str(sequences[seq2].letter_annotations["phred_quality"])
    pairCount = pairCount + 1
    #print(alignment)
    #print("\n ******************************** \n")
    print(alignment[pairCount])
    print("\n ******************************** \n")
"""
# Scoring function
for key in alignment:
    alignment[key]
    
def overalpScore(seq1, seq2, i1, i2, L):
   """ 







###########################################################################
"""
 #TRIALS
 
 
#print(pairwise2.format_alignment(*al[0]))
#print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n")
al = pairwise2.align.globalxx("AGTCGCTCGTACGCTACCGAT","AGTCGCTCGCGTATACCGAT")
#print(al[0])
#a = pairwise2.format_alignment(*al)
print(al)
print("\n +")
print(al[0][1])
#result.write(pairwise2.format_alignment(al)) 
#
#print(al)

#alignments = pairwise2.align.globalxx("ACCGT", "ACG")
#print(alignments[1])
#print(pairwise2.format_alignment(*alignments))

"""
