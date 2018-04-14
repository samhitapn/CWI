#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 14:26:16 2018
@author: samhitapn
@purpose: Overlap score trials
@modified: 12.April.2018
"""

#import biopython
from Bio import pairwise2 
from Bio import SeqIO as sio
import numpy
import itertools
from time import time

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


"""
@Definition: For calculating the probaility from the phred quality scores
@Input parameters: 
@Output parameters: Probability
"""
def getProbQuality (q):
    #print(q)
    q = int(q)
    p = 10**(-q/10)
    return p


nt = ["A","T","C","G"] # Bases
"""
@Definition: 
@Input parameters: 
@Output parameters: 
"""   
#def overalpScoreCalculation(seqDetails, i1, i2, L):
def overalpScoreCalculation(seqDetails, i):
    probabilityOverall = 1
    scoresRead1 = seqDetails.split(";")[2].split(",")
    scoresRead2 = seqDetails.split(";")[5].split(",")
    
    L = abs(len(seqDetails.split(";")[4]) - i)
    print(len(seqDetails.split(";")[4]),len(seqDetails.split(";")[2]),i,L)
    i1 = i
    i2 = 0
    while i1 < i1 + L and i2 < i2 + L:
        probabilityBase = 0
        
        #New code -> To include indels Option1
        #scoresRead1 = seqDetails.split(";")[2].split(",")
        #scoresRead2 = seqDetails.split(";")[5].split(",")
            # Gap in first read -> calculation based on read 2
        if seqDetails.split(";")[1][i1] == "-":
           #print(scoresRead2[i2])
           probabilityBase = (3/13 * getProbQuality(float(scoresRead2[i2]))) + (10/13 * (1 - getProbQuality(float(scoresRead2[i2]))))
           pl = 1
           # Gap in second read -> calculation based on read 1
        elif seqDetails.split(";")[4][i2] == "-":
           print(scoresRead1[i1])
           probabilityBase = (3/13 * getProbQuality(float(scoresRead1[i1]))) + (10/13 * (1 - getProbQuality(float(scoresRead1[i1]))))
           pl = 2
           # Existing score calculation
        else:
            for n in nt:
                #print(n)
                probabilityBase = probabilityBase + (probabilityQ(n,seqDetails.split(";")[1][i1],getProbQuality(float(scoresRead1[i1]))) * probabilityQ(n,seqDetails.split(";")[4][i2],getProbQuality(float(scoresRead2[i2]))))
            pl = 3
        probabilityOverall = probabilityOverall * probabilityBase
        print(i1,i2,probabilityOverall,pl)
        i1 = i1 + 1
        i2 = i2 + 1
    # Overlap score
    overlapScore = probabilityOverall ** 1/L
    #print(overlapScore)
    return (overlapScore)
 
# MAIN


#overlapFile = f.readlines()
#print(type(overlapFile))
#print("^^^^^OVERLAPS READ^^^^^^^\n")
# Read in the fastq files -> Usually the overalp pairs; here only the test 5 reads
sequences = sio.to_dict(sio.parse("data/Sample_AllReads.fastq","fastq"))
print("^^^^^^SEQUENCES READ^^^^^^\n")

overlapCount = 1
alignment = {}
# Get the overlap pairs and details from the PAF files
with open("data/Sample_AllReads_Overlaps.paf","r") as f:
    for ovlPair in f.readlines():
        if overlapCount <= 30:
        #print(ovlPair)
        #print("^^^^^^^^^^^^^^^^^^^^")
            
            ovl = ovlPair.split("\t")
            if ovl[0] != ovl[5]:
                print(overlapCount)
                start = time()
                al = pairwise2.align.globalxx(sequences[ovl[0]].seq,sequences[ovl[5]].seq)
                stop = time()
                print(stop - start)
                alignment[ovl[0] + ' ' + ovl[5]] = ovl[10] + " " + al[0][0] + " " + str(sequences[ovl[0]].letter_annotations["phred_quality"]).strip('[]').replace(" ","") + " " + al[0][1] + " " + str(sequences[ovl[0]].letter_annotations["phred_quality"]).strip('[]').replace(" ","")
                overlapCount = overlapCount + 1
        else:
            break
print(alignment)


"""
# Read in the fastq files -> Usually the overalp pairs; here only the test 5 reads
sequences = sio.to_dict(sio.parse("data/Sample_AllReads.fastq","fastq"))
for key, value in sequences.items():
    print(key)
    print("\n ^^^^ \n")


# Get the alignment of all the overlap pairs
alignment = {}
pairCount = 1
for seq1, seq2 in itertools.combinations(sequences, 2):
    al = pairwise2.align.globalxx(sequences[seq1].seq,sequences[seq2].seq)
    alignment[sequences[seq1].id + '&' + sequences[seq2].id] = sequences[seq1].id + ";" + al[0][0] + ";" + str(sequences[seq1].letter_annotations["phred_quality"]).strip('[]').replace(" ","") + ";" + sequences[seq2].id + ";" + al[0][1] +  ";" + str(sequences[seq2].letter_annotations["phred_quality"]).strip('[]').replace(" ","")
    #print(alignment)
    #print(alignment[pairCount].split(",")[1][3])
    #print("\n ******************************** \n")
    #pairCount = pairCount + 1


# Getting the scores for each overlap pair
for key in alignment:
    print(alignment[key].split(";")[1])
    print(alignment[key].split(";")[4])
    overlapScore = overalpScoreCalculation(alignment[key],20)
    
    alignment[key] = alignment[key] + ";" + str(overlapScore)
    print(alignment[key].split(";")[6])
    print("\n ############ \n")
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
