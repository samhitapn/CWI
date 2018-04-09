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
    #print(alignment)
    #print(alignment[pairCount].split(",")[1][3])
    #print("\n ******************************** \n")
    pairCount = pairCount + 1


"""
@Definition: For calculating the probaility of the base at a given position
@Input parameters: X -> base A/C/T/G; b -> base in the particular position; q -> quality score at the position
@Output parameters: Probability
"""
def probabilityQ (X, b, p):
    if X == b:
        qp = 1 - p
    else:
        qp = p/3
    return qp



nt = ["A","T","C","G"] # Bases

"""
@Definition: For calculating the probaility of the base at a given position
@Input parameters: X -> base A/C/T/G; b -> base in the particular position; q -> quality score at the position
@Output parameters: Probability
"""   
def overalpScoreCalculation(seqDetails, i1, i2, L):
    probabilityOverall = 1
    while i1 <= i1 + L and i2 <= i2 + L:
        probabilityBase = 0
        
        #New code -> To include indels Option1
        
            # Gap in first read -> calculation based on read 2
        print(seqDetails)
        if seqDetails.split(",")[1][i1] == "-":
           probabilityBase = (3/13 * float(seqDetails.split(",")[5][i2])) + (10/13 * (1 - float(seqDetails.split(",")[5][i2])))
           # Gap in second read -> calculation based on read 1
        elif seqDetails.split(",")[4][i2] == "-":
           probabilityBase = (3/13 * float(seqDetails.split(",")[2][i1])) + (10/13 * (1 - float(seqDetails.split(",")[2][i1])))
           # Existing score calculation
        else:
            for n in nt:
                probabilityBase = probabilityBase + (probabilityQ(n,seqDetails.split(",")[1][i1],float(seqDetails.split(",")[2][i1])) * probabilityQ(n,seqDetails.split(",")[4][i2],float(seqDetails.split(",")[5][i2])))
        probabilityOverall = probabilityOverall * probabilityBase
        i1 = i1 + 1
        i2 = i2 + 1
    # Overlap score
    overlapScore = probabilityOverall ** 1/L
    print(overlapScore)
    return (overlapScore)
 
# Scoring function
for key in alignment:
    #print(alignment[key])
    overlapScore = overalpScoreCalculation(alignment[key],20,1,10)
    alignment[key] = alignment[key] + "," + overlapScore
    print(alignment[key].split(",")[0],alignment[key].split(",")[3],alignment[key].split(",")[6])
    print("\n ############ \n")
    

    
        



     



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
