#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: 02.may.2018
@author: samhitapn
@purpose: Overlap Score generator along with preprocessing
@modified: 02.May.2018
"""

# LIBRARIES IMPORT
from Bio import pairwise2
from Bio import SeqIO as sio
import numpy as np
import math as mt
import itertools
from time import time
import re
import os
import sys
from Bio.Seq import MutableSeq

# GLOBAL VARIABLES REQUIRED
cigarPattern = re.compile('([0-9]*)([IDM])')
nt = ["A","T","C","G"] # Bases

# FUNCTIONS
"""
@Definition: Mapping the sequences, getting the sequence and the scores and converting the CIGAR string to the alignment
@Input parameters:
@Output parameters:
"""
def getSeqFromCigar(cigar):
    cigarSeq = []
    for num, char in cigarPattern.findall(cigar):
        if num:
            num = int(num)
        else:
            num = 1
        cigarSeq.append("".join(char * num))
    cigarSeq = "".join(cigarSeq)
    return(cigarSeq)

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
    #q = int(q)
    p = 10**(-q/10)
    return p

"""
@Definition:
@Input parameters:
@Output parameters: Probability
"""
def getGapRegionScore (score, lenGap):
    probSum = 0
    #if lenGap > 1:
    for i in range(0, len(score)):
            probSum = probSum + (10/13 * getProbQuality(score[i])) + (3/13 * (1 - getProbQuality(score[i])))
    probGapSum = probSum/lenGap
    return(probGapSum)

"""
@Definition: Mapping the sequences, getting the sequence and the scores and converting the CIGAR string to the alignment
@Input parameters:
@Output parameters:
"""
#from signal import signal, SIGPIPE, SIG_DFL
#signal(SIGPIPE, SIG_DFL)
def getOverlapScore(readData):
    probabilityBase = 0
    probabilityOverall = 1
    cigarSeq = getSeqFromCigar(readData[8])
    print(len(cigarSeq))
    seq1 = readData[0][readData[2]:]
    seq2 = readData[4][readData[6]:]
    score1 = readData[1][readData[2]:]
    score2 = readData[5][readData[6]:]
    print(readData[2],readData[3],readData[6],readData[7])
    print(len(seq1),len(seq2),len(score1),len(score2))
    pos1 = 0
    pos2 = 0
    L = 0
    for num, char in cigarPattern.findall(readData[8]):
        #pos = 0
        #print(num,char)

        if num:
            num = int(num)
        else:
            num = 1
        L = L + num
        if char == "I":
            probabilityBase = probabilityBase + getGapRegionScore(score2[pos2:pos2+num], num)
            pos2 = pos2 + num
            #print(pos,char,probabilityBase)
        if char == "D":
            probabilityBase = probabilityBase + getGapRegionScore(score1[pos1:pos1+num], num)
            pos1 = pos1 + num
            #print(pos,char,probabilityBase)
        if char == "M":
            tempSeq1 = seq1[pos1:pos1 + num]
            tempScore1 = score1[pos1:pos1 + num]
            tempSeq2 = seq2[pos2:pos2 + num]
            tempScore2 = score2[pos2:pos2 + num]
            for i in range(0, len(tempSeq1)):
                for n in nt:
                    probabilityBase = probabilityBase + (probabilityQ(n,tempSeq1[i],getProbQuality(tempScore1[i])) * probabilityQ(n,tempSeq2[i],getProbQuality(tempScore2[i])))
            pos1 = pos1 + num
            pos2 = pos2 + num
        #print(num,pos,char,probabilityBase)
        probabilityOverall = probabilityOverall * probabilityBase
    overlapScore = probabilityOverall ** 1/L
    return(overlapScore)



"""
    for i in range(0,len(cigarSeq)):
        if cigarSeq[i] == "I":
            #seq1.insert(i,"-")
            #score1.insert(i,"-")

        elif cigarSeq[i] == "D":
            seq2.insert(i,"-")
            score2.insert(i,"-")
    readData[0] = seq1
    readData[1] = score1
    readData[4] = seq2
    readData[5] = score2
    return(readData)
"""

# MAIN
readPairData = dict()
with open("data/sequences/fastq_100Reads_NewNames/100Reads_All_Overlaps.paf") as paf:
    pafData = paf.readlines()

fastq = sio.to_dict(sio.parse("data/sequences/fastq_100Reads_NewNames/fastq_merged_100Reads.fastq","fastq"))

for overlapPair in pafData:
    tempData = list()
    ovl = overlapPair.split("\t")
    print(overlapPair)
    tempData = [list(fastq[ovl[0]].seq), fastq[ovl[0]].letter_annotations["phred_quality"], int(ovl[2]), int(ovl[3]), list(fastq[ovl[5]].seq),fastq[ovl[5]].letter_annotations["phred_quality"],int(ovl[7]), int(ovl[8]),ovl[20].split(":")[2].strip()]
    readPairData[ovl[0] + "-" + ovl[5]] = tempData
    break
#print(readPairData)

for key in readPairData:
    print(getOverlapScore(readPairData[key]))
    break
