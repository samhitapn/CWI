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
    #newSeq1 = readData[0]
    #newSeq2 = readData[4]
    #newScore1 = readData[1]
    #newScore2 = readData[5]
    start1 = readData[2]
    start2 = readData[6]
    for i in range(0,len(cigarSeq)):
        print(cigarSeq[i])
        if cigarSeq[i] == "I":
            readData[0][start1] = "-"
            readData[1].insert(start1,"-")
        elif cigarSeq[i] == "D":
            readData[4][start2] = "-"
            readData[5].insert(start2,"-")
        start1 = start1 + 1
        start2 = start2 + 1
    return(readData)


# MAIN
readPairData = dict()
with open("data/sequences/fastq_100Reads_NewNames/100Reads_All_Overlaps.paf") as paf:
    pafData = paf.readlines()

fastq = sio.to_dict(sio.parse("data/sequences/fastq_100Reads_NewNames/fastq_merged_100Reads.fastq","fastq"))

for overlapPair in pafData:
    tempData = list()
    ovl = overlapPair.split("\t")
    #print(fastq[ovl[0]].seq)
    tempData = [list(fastq[ovl[0]].seq), str(fastq[ovl[0]].letter_annotations["phred_quality"]).strip('[]').replace(" ",""), int(ovl[2]), int(ovl[3]), list(fastq[ovl[5]].seq), str(fastq[ovl[5]].letter_annotations["phred_quality"]).strip('[]').replace(" ",""), int(ovl[7]), int(ovl[8]),ovl[20]]
    readPairData[ovl[0] + "-" + ovl[5]] = tempData
#print(readPairData)

for key in readPairData:
    #print(key)
    #print(readPairData[key][0])
    #print("&&&&&&&&")
    #break
    finalAlignment = getFinalAlignment(readPairData[key])
