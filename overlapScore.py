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
@Definition: Mapping the sequences, getting the sequence and the scores and converting the CIGAR string to the alignment
@Input parameters:
@Output parameters:
"""
#from signal import signal, SIGPIPE, SIG_DFL
#signal(SIGPIPE, SIG_DFL)
def getFinalAlignment(readData):
    cigarSeq = getSeqFromCigar(readData[8])
    seq1 = readData[0][readData[2]:readData[3]]
    seq2 = readData[4][readData[6]:readData[7]]
    score1 = readData[1][readData[2]:readData[3]]
    score2 = readData[5][readData[6]:readData[7]]

    for i in range(0,len(cigarSeq)):
        if cigarSeq[i] == "I":
            seq1.insert(i,"-")
            score1.insert(i,"-")
        elif cigarSeq[i] == "D":
            seq2.insert(i,"-")
            score2.insert(i,"-")
    readData[0] = seq1
    readData[1] = score1
    readData[4] = seq2
    readData[5] = score2
    return(readData)

"""
@Definition:
@Input parameters:
@Output parameters:
"""
def overalpScoreCalculation(seqDetails):
    # Initiating required values
    probabilityOverallSum = 1

    # Getting the sequence and scores
    seqRead1 =  seqDetails[0]
    scoreRead1 = seqDetails[1]
    seqRead2 = seqDetails[4]
    scoreRead2 = seqDetails[5]

    # Gap in first read -> calculation based on read 2



# MAIN
readPairData = dict()
with open("data/sequences/fastq_100Reads_NewNames/100Reads_All_Overlaps.paf") as paf:
    pafData = paf.readlines()

fastq = sio.to_dict(sio.parse("data/sequences/fastq_100Reads_NewNames/fastq_merged_100Reads.fastq","fastq"))

for overlapPair in pafData:
    tempData = list()
    ovl = overlapPair.split("\t")
    #print(fastq[ovl[0]].seq)
    tempData = [list(fastq[ovl[0]].seq), fastq[ovl[0]].letter_annotations["phred_quality"], int(ovl[2]), int(ovl[3]), list(fastq[ovl[5]].seq),fastq[ovl[5]].letter_annotations["phred_quality"],int(ovl[7]), int(ovl[8]),ovl[20].split(":")[2].strip()]
    readPairData[ovl[0] + "-" + ovl[5]] = tempData
#print(readPairData)

for key in readPairData:
    #print(len(readPairData[key][0]),len(readPairData[key][1]),len(readPairData[key][4]),len(readPairData[key][5]))
    readPairData[key] = getFinalAlignment(readPairData[key])
    #print(len(readPairData[key][0]),len(readPairData[key][1]),len(readPairData[key][4]),len(readPairData[key][5]))

    print(readPairData[key][0],readPairData[key][4])
    #print("^^^^^^^^^^^")
