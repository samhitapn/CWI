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
    print(s.letter_annotations["phred_quality"][1])
    #print(s.format("fastq"))
"""
sequences = sio.to_dict(sio.parse("data/ONT_Sample1_5Reads.fastq","fastq"))
#print(seq)


alignment = {}
pairCount = 1
for seq1, seq2 in itertools.combinations(sequences, 2):
    al = pairwise2.align.globalxx(sequences[seq1].seq,sequences[seq2].seq)
    alignment[pairCount] = sequences[seq1].id + "," + sequences[seq2].id + "," + al[0][0] + "," + str(sequences[seq1].letter_annotations["phred_quality"]) + "," + al[0][1] +  "," + str(sequences[seq2].letter_annotations["phred_quality"])
    pairCount = pairCount + 1
    print(alignment)
    print("\n ******************************** \n")
    

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
