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
def getFinalAlignment(readData):
    cigarSeq = getSeqFromCigar(readData[8])
    newSeq1 = []
    newSeq2 = []
    for i in range(readData[2],readData[3] + 1):
        print(i,readData[2][i])


# MAIN
readPairData = dict()
with open("data/sequences/fastq_100Reads_NewNames/100Reads_All_Overlaps.paf") as paf:
    pafData = paf.readlines()

fastq = sio.to_dict(sio.parse("data/sequences/fastq_100Reads_NewNames/fastq_merged_100Reads.fastq","fastq"))

for overlapPair in pafData:
    tempData = list()
    ovl = overlapPair.split("\t")
    #print(fastq[ovl[0]].seq)
    tempData = [fastq[ovl[0]].seq, str(fastq[ovl[0]].letter_annotations["phred_quality"]), int(ovl[2]), int(ovl[3]), fastq[ovl[5]].seq, str(fastq[ovl[5]].letter_annotations["phred_quality"]), int(ovl[7]), int(ovl[8]),ovl[20]]
    readPairData[ovl[0] + "-" + ovl[5]] = tempData
#print(readPairData)

for key in readPairData:
    #print(key)
    #print(readPairData[key][0])
    #print("&&&&&&&&")
    #break
    finalAlignment = getFinalAlignment(readPairData[key])
