#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: 02.may.2018
@author: samhitapn
@purpose: Overlap Score generator along with preprocessing
@modified: 02.May.2018
"""

# LIBRARIES IMPORT
#from Bio import pairwise2
#from Bio import SeqIO as sio
import numpy as np
import math as mt
import itertools
from time import time
import re
import os
import sys
import argparse


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
        #print("Same",qp)
    else:
        qp = pNew/3
        #print("diff",qp)
    return qp


"""
@Definition: For calculating the probaility from the phred quality scores
@Input parameters:
@Output parameters: Probability
"""
def getProbQuality (q):
    #print(q)
    #q = int(q)
    qNew = ord(q)  - 33
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
def getGapRegionScore (score, lenGap, scoreType):
    gapProb = 0
    #print(scoreList)
    #scoreList = []
    #if lenGap > 1:
    convertedScores = [getProbQuality(score[i]) for i in range(0,lenGap)]
    score = lambda sc : (97/100 * sc) + (3/100 * (1 - sc))
    """
    for i in range(0, lenGap):
        #print(score[i])
        #score[i] = getProbQuality(score[i])
        #scoreList = scoreList.append(getProbQuality(score[i]))
        #print("Score&&&&&",ord(score[i]))
        ovlScore = ovlScore + getProbQuality(score[i])
    ovlScore = ovlScore/lenGap
    """
    if scoreType == 0:
        gapScore = max(convertedScores)
        gapProb = score(gapScore)
    elif scoreType == 1:
        gapScore = min(convertedScores)
        gapProb = score(gapScore)
    elif scoreType == 2:
        gapScore = np.mean(convertedScores)
        gapProb = score(gapScore)
    elif scoreType == 3:
        gapScore = np.median(convertedScores)
        gapProb = score(gapScore)
    elif scoreType == 4:
        gapScore = np.prod(convertedScores) ** (1/len(convertedScores))
        gapProb = score(gapScore)
    elif scoreType == 5:
        for i in convertedScores:
            gapProb = gapProb + score(i)
    #gapProb = (10/13 * ovlScore) + (3/13 * (1 - ovlScore))
    #gapProb = (97/100 * gapScore) + (3/100 * (1 - gapScore))
    #print("PROB", prob)
    #print(scoreList)
    return(gapProb)

"""
@Definition: Mapping the sequences, getting the sequence and the scores and converting the CIGAR string to the alignment
@Input parameters:
@Output parameters:
"""

def getOverlapScore(key, readData, scoreType):
    #scoreList = [0]
    probabilityOverall = 0
    probabilityGaps = 0
    probabilityMatches = 0
    prob = 0
    seq1 = readData[0]
    seq2 = readData[5]
    score1 = readData[1]
    score2 = readData[6]
    result = []
    alpha1 = (readData[4] - readData[3])/(readData[2] - readData[3] + readData[8])
    alpha2 = (readData[9] - readData[8])/(readData[9] + readData[2] - readData[4])
    alpha = max(alpha1, alpha2)
    #c = 0
    #cig = readData[8].count("M") + readData[8].count("I") + readData[8].count("D")
    #print("CIGARDetails:",cig, readData[8].count("M"), readData[8].count("I"), readData[8].count("D"))
    #print("CIGAR:",readData[8])
    #print(len(seq1),len(seq2))
    for i in range(0,3):
        print(i)
        pos1 = readData[3]
        pos2 = readData[8]
        L = 0
        for num, char in cigarPattern.findall(readData[10]):

            if num:
                num = int(num)
            else:
                num = 1
            #c = c + num
            #if pos1 <= readData[3] and pos2 <= readData[7]:
            tempSeq1 = seq1[pos1:pos1 + num]
            tempScore1 = score1[pos1:pos1 + num]
            tempSeq2 = seq2[pos2:pos2 + num]
            tempScore2 = score2[pos2:pos2 + num]

            if char == "I":
                if i == 1:
                    pos1 = pos1 + num
                else:
                    sc = getGapRegionScore(tempScore1, num, scoreType)
                    prob = prob + np.log(sc)
                    pos1 = pos1 + num
                    L = L + 1
            elif char == "D":
                if i == 1:
                    pos2 = pos2 + num
                else:
                    sc = getGapRegionScore(tempScore2, num, scoreType)
                    prob = prob + np.log(sc)
                    pos2 = pos2 + num
                    L = L + 1
            elif char == "M":
                if i == 0:
                    pos1 = pos1 + num
                    pos2 = pos2 + num
                else:
                    for i in range(0, num):
                        probabilityBase = 0
                        for n in nt:
                            sc = (probabilityQ(n,tempSeq1[i],tempScore1[i]) * probabilityQ(n,tempSeq2[i],tempScore2[i]))
                            probabilityBase = probabilityBase + sc
                            #probabilityBase = probabilityBase + (probabilityQ(n,tempSeq1[i],np.float128(tempScore1[i])) * probabilityQ(n,tempSeq2[i],np.float128(tempScore2[i])))
                        prob = prob + np.log(probabilityBase)
                        L = L + 1
                    pos1 = pos1 + num
                    pos2 = pos2 + num
                    #print("M:",num,pos1,pos2)
        print(alpha1,alpha2,alpha)
        overlapScore = (np.exp(prob) ** (1/L))
        result.append(overlapScore)
        print(result)
    return(result)


# MAIN
# Arguement parsing
parser = argparse.ArgumentParser(description='Parser for input files and output files')
parser.add_argument('-p','--paf', help='PAF file name',required=True)
parser.add_argument('-f','--fastq',help='Fastq file name', required=True)
parser.add_argument('-o','--output',help='Fastq file name', required=True)
parser.add_argument('-g','--gapScoreType',help='Gap region method max (0), min (1), avg (2), median (3), geometricMean (4), position-by-position (5), all(6)', required=True, type = int)
#parser.add_argument('-s','--scoreRegion',help='Score only gaps and/or substitutions: gaps (0), substitutions (1), both (2), all (3)', required=True, type = int)
args = parser.parse_args()

readPairData = dict()
fastqTemp = {}
with open(args.paf) as paf:
#with open("data/sequences/fastq_100Reads_NewNames/sample.paf") as paf:
    pafData = paf.readlines()

#fastq = sio.to_dict(sio.parse("data/sequences/fastq_100Reads_NewNames/fastq_merged_100Reads.fastq","fastq"))
with open(args.fastq) as fastq:
#with open("data/sequences/fastq_100Reads_NewNames/sample.fastq") as fastq:
    fastqData = fastq.readlines()

for i in range(0,len(fastqData)):
    if fastqData[i].startswith("@"):
         fastqTemp[fastqData[i].split(" ")[0].strip("@")] = [fastqData[i+1].rstrip("\n"),fastqData[i+3].rstrip("\n")]
    i = i + 1

statFileName = args.output + "_stats.csv"
statFile = open(statFileName,"w+")
statFile.write("KEY \t READ1_LENGTH \t READ1_START \t READ1_END \t READ1_OVL_LEN \t READ2_LENGTH \t READ2_START \t READ2_END \t READ2_OVL_LEN \t MATCHES \t INSERTIONS \t DELETIONS \n")
for overlapPair in pafData:
    tempData = list()
    ovl = overlapPair.split("\t")

    """
    if ovl[0] + "-" + ovl[5] == "origSeq_159-origSeq_77":

        print(ovl)
        print(fastqTemp[ovl[0]][0])
        print(fastqTemp[ovl[0]][1])
        print(fastqTemp[ovl[5]][0])
        print(fastqTemp[ovl[5]][1])

        print(ovl)
        print("Init:",ovl[20].split(":")[2].strip("\n"))
    """
    cigarIndex = [ovl.index(i) for i in ovl if i.startswith("cg")]
    cig = getSeqFromCigar(ovl[cigarIndex[0]].split(":")[2].strip("\n"))
    statFile.write(str(ovl[0] + "-" + ovl[5]) + "\t" + str(ovl[1]) + "\t" + str(ovl[2]) + "\t" + str(ovl[3]) + "\t" + str(int(ovl[3]) - int(ovl[2])) + "\t" + str(ovl[6]) + "\t" + str(ovl[7]) + "\t" + str(ovl[8]) + "\t" + str(int(ovl[8]) - int(ovl[7])) + "\t" + str(cig.count("M")) + "\t" + str(cig.count("I")) + "\t" + str(cig.count("D")) + "\n")
    readPairData[ovl[0] + "-" + ovl[5]] = [fastqTemp[ovl[0]][0],fastqTemp[ovl[0]][1], int(ovl[1]), int(ovl[2]), int(ovl[3]), fastqTemp[ovl[5]][0], fastqTemp[ovl[5]][1], int(ovl[6]), int(ovl[7]),int(ovl[8]),ovl[cigarIndex[0]].split(":")[2].strip("\n")]
    #tempData = [list(fastq[ovl[0]].seq), fastq[ovl[0]].letter_annotations["phred_quality"], int(ovl[2]), int(ovl[3]), list(fastq[ovl[5]].seq),fastq[ovl[5]].letter_annotations["phred_quality"],int(ovl[7]), int(ovl[8]),ovl[20].split(":")[2].strip()]
    #readPairData[ovl[0] + "-" + ovl[5]] = [fastqTemp[ovl]]
statFile.close()

print(len(readPairData))
c = 0
#scoreFileName = args.output + "_scores.csv"
#outputFile = open(scoreFileName,"w+")
"""
for key in readPairData:
    keyElements = key.split("-")
    if keyElements[0].split("_")[0] == keyElements[1].split("_")[0]:
        ovlType = "Good Overlap"
    else:
        ovlType = "Bad Overlap"

    if args.gapScoreType == 6:
        gapScoreType = [0,1,2,3,4,5]
    else:
        gapScoreType = [args.gapScoreType]
"""
if args.gapScoreType == 6:
    gapScoreType = [0,1,2,3,4,5]
else:
    gapScoreType = [args.gapScoreType]

for i in gapScoreType:
    if i == 0:
        name = "max"
        print(name)
    elif i == 1:
        name = "min"
        print(name)
    elif i == 2:
        name = "mean"
        print(name)
    elif i == 3:
        name = "median"
        print(name)
    elif i == 4:
        name = "geometricMean"
        print(name)
    elif i == 5:
        name = "positionByPosition"
        print(name)

    scoreFileName = args.output + name + "_scores.csv"
    outputFile = open(scoreFileName,"w+")
    outputFile.write("KEY \t GAPS \t MATCHES \t ALL \t OVERLAP_TYPE \n")

    for key in readPairData:
        keyElements = key.split("-")
        if keyElements[0].split("_")[0] == keyElements[1].split("_")[0]:
            ovlType = "Good Overlap"
        else:
            ovlType = "Bad Overlap"

        results = getOverlapScore(key, readPairData[key],i)

        """
        if results[0] == 1:
            error = "IndexError"
        else:
            error = "No error"
        """
        print(c,key,str(results[0]))

        outputFile.write(key + "\t" + str(results[0]) + "\t" + str(results[1]) + "\t" + str(results[2]) + "\t" + ovlType + "\n")
    outputFile.close()
