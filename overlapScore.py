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
    pNew = getProbQuality(p)
    if X == b:
        qp = 1 - pNew
    else:
        qp = pNew/3
    return qp


"""
@Definition: For calculating the probaility from the phred quality scores
@Input parameters:
@Output parameters: Probability
"""
def getProbQuality (q):
    #print(q)
    #q = int(q)
    qNew = ord(q)
    #print(q,qNew)
    #p = 10**(-np.float128(q)/10)
    p = 10 ** (-np.float128(qNew)/10)
    #print(p)
    return p

"""
@Definition:
@Input parameters:
@Output parameters: Probability
"""
def getGapRegionScore (score, lenGap):
    ovlScore = 0
    #if lenGap > 1:

    for i in range(0, lenGap):
        #print(score[i])
        #score[i] = getProbQuality(score[i])
        ovlScore = ovlScore + getProbQuality(score[i])
    ovlScore = ovlScore/lenGap
    prob = (10/13 * ovlScore) + (3/13 * (1 - ovlScore))
    return(prob)

"""
@Definition: Mapping the sequences, getting the sequence and the scores and converting the CIGAR string to the alignment
@Input parameters:
@Output parameters:
"""
#from signal import signal, SIGPIPE, SIG_DFL
#signal(SIGPIPE, SIG_DFL)
def getOverlapScore(key, readData):
    #probabilityBase = 0
    probabilityOverall = 1
    #cigarSeq = getSeqFromCigar(readData[8])
    #print(len(cigarSeq))
    #seq1 = readData[0][readData[2]:]
    #seq2 = readData[4][readData[6]:]
    #score1 = readData[1][readData[2]:]
    #score2 = readData[5][readData[6]:]
    seq1 = readData[0]
    seq2 = readData[4]
    score1 = readData[1]
    score2 = readData[5]
    #print(readData[2],readData[3],readData[6],readData[7])
    #print(len(seq1),len(seq2),len(score1),len(score2))
    pos1 = readData[2]
    pos2 = readData[6]
    L = 0
    errorDet = list()
    for num, char in cigarPattern.findall(readData[8]):
            #pos = 0
           #print(num,char)
        try:
            if num:
                num = int(num)
            else:
                num = 1
            #L = L + num
            #print(key, num, char, pos1, pos1 + num, pos2, pos2 + num,len(seq1),len(seq2))
            tempSeq1 = seq1[pos1:pos1 + num]
            tempScore1 = score1[pos1:pos1 + num]
            tempSeq2 = seq2[pos2:pos2 + num]
            tempScore2 = score2[pos2:pos2 + num]

            if char == "I":
                #probabilityBase = probabilityBase + getGapRegionScore(score2[pos2:pos2+num], num)
                #print("I",probabilityOverall)
                probabilityOverall = probabilityOverall * getGapRegionScore(tempScore2, num)
                print(char,num,probabilityOverall)
                assert 0 <= probabilityOverall <= 1, print(char, pos1,pos2,probabilityOverall)
                pos2 = pos2 + num
                L = L + 1
                #print(num,pos1,pos2,char,probabilityBase)
            elif char == "D":
                #probabilityBase = probabilityBase + getGapRegionScore(score1[pos1:pos1+num], num)
                probabilityOverall = probabilityOverall * getGapRegionScore(tempScore1, num)
                print(char,num,probabilityOverall)
                assert 0 <= probabilityOverall <= 1, print(char, pos1,pos2,probabilityOverall)
                pos1 = pos1 + num
                L = L + 1
                #print(num,pos1,pos2,char,probabilityBase)
            elif char == "M":
                for i in range(0, num):
                    probabilityBase = 0
                    print(tempScore1[i],getProbQuality(tempScore1[i]),tempScore2[i],getProbQuality(tempScore2[i]))
                    for n in nt:
                        probabilityBase = probabilityBase + (probabilityQ(n,tempSeq1[i],tempScore1[i]) * probabilityQ(n,tempSeq2[i],tempScore2[i]))
                        assert 0 <= probabilityBase <= 1, print(char, pos1,pos2,probabilityBase, "Base")
                        #probabilityBase = probabilityBase + (probabilityQ(n,tempSeq1[i],np.float128(tempScore1[i])) * probabilityQ(n,tempSeq2[i],np.float128(tempScore2[i])))
                    #assert 0 < probabilityOverall < 1, print(char, pos1,pos2,probabilityOverall)
                    probabilityOverall = probabilityOverall * probabilityBase
                    print(char,num,probabilityBase,probabilityOverall)
                    assert 0 <= probabilityOverall <= 1, print(char, pos1,pos2,probabilityOverall)
                    L = L + 1
                pos1 = pos1 + num
                pos2 = pos2 + num
                #print(num,pos1,pos2,char,probabilityBase)
            #probabilityOverall = probabilityOverall * probabilityBase
            #print(probabilityOverall)
        except IndexError:
            result = [pos1, pos2, 0]
            continue
        else:
            overlapScore = probabilityOverall ** (1/L)
            print(L, probabilityOverall,overlapScore)
            result = [pos1, pos2, overlapScore]
        #break
    return(result)


# MAIN
readPairData = dict()
fastqTemp = {}
with open("data/sequences/fastq_100Reads_NewNames/100Reads_All_Overlaps.paf") as paf:
    pafData = paf.readlines()

#fastq = sio.to_dict(sio.parse("data/sequences/fastq_100Reads_NewNames/fastq_merged_100Reads.fastq","fastq"))
with open("data/sequences/fastq_100Reads_NewNames/fastq_merged_100Reads.fastq") as fastq:
    fastqData = fastq.readlines()

for i in range(0,len(fastqData)):
    if fastqData[i].startswith("@"):
         fastqTemp[fastqData[i].split(" ")[0].strip("@")] = [fastqData[i+1].rstrip("\n"),fastqData[i+3].rstrip("\n")]
    i = i + 1

for overlapPair in pafData:
    tempData = list()
    ovl = overlapPair.split("\t")
    #print(overlapPair)
    """
    if str(ovl[0] + "-" + ovl[5])  == "seq5_16-seq8_86":
        print(overlapPair)
        print("".join(list(fastq[ovl[0]].seq)))
        print("".join(str(fastq[ovl[0]].letter_annotations["phred_quality"])))
        print("****")
        print("".join(list(fastq[ovl[5]].seq)))
        print("".join(str(fastq[ovl[5]].letter_annotations["phred_quality"])))
    """
    readPairData[ovl[0] + "-" + ovl[5]] = [fastqTemp[ovl[0]][0],fastqTemp[ovl[0]][1],int(ovl[2]), int(ovl[3]),fastqTemp[ovl[5]][0],fastqTemp[ovl[5]][1],int(ovl[7]),int(ovl[8]),ovl[20].split(":")[2].strip("\n")]
    #tempData = [list(fastq[ovl[0]].seq), fastq[ovl[0]].letter_annotations["phred_quality"], int(ovl[2]), int(ovl[3]), list(fastq[ovl[5]].seq),fastq[ovl[5]].letter_annotations["phred_quality"],int(ovl[7]), int(ovl[8]),ovl[20].split(":")[2].strip()]
    #readPairData[ovl[0] + "-" + ovl[5]] = [fastqTemp[ovl]]

print(len(readPairData))
#print(readPairData)
c = 0
outputFile = open("100ReadsOutput.txt","w+")
test = open("test.txt","w+")
for key in readPairData:
    #if key == "seq5_16-seq8_86":
    print(readPairData[key])
    for i in readPairData[key]:
        test.write(readPairData[key][i])
    c = c + 1
    results = getOverlapScore(key, readPairData[key])
    #print(results)
    keyElements = key.split("-")
    if keyElements[0].split("_")[0] == keyElements[1].split("_")[0]:
        ovlType = "Good Overlap"
    else:
        ovlType = "Bad Overlap"

    if results[2] == 0:
        error = "IndexError"
    else:
        error = "No error"
    #readPairData[key].append(score)
    print(key,str(results[2]))
    outputFile.write(key + "\t" + error + "\t" + str(results[2]) + "\t" + ovlType + "\n")
    break
outputFile.close()
test.close()

    #print(c, key, readPairData[key][9])
