#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: 14.aug.2018
@author: samhitapn
@purpose: Overlap Score generator along with preprocessing ignoring the presence of gaps
@modified: 14.Aug.2018
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
@Definition: Converting the CIGAR string to the expanded form
@Input parameters: CIGAR string
@Output parameters: expanded CIGAR string
"""
def getExpandedCigar(cigar):
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
@Definition: For calculating the probability of the base at a given position
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
@Input parameters: Raw phred scores
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
@Definition: For calculating the probability of the gap region
@Input parameters: List of scores, length of the gap region, scoring type
@Output parameters: Probability
"""
def getGapRegionScore (score, lenGap, scoreType):
    gapProb = 0
    #print(scoreList)
    #scoreList = []
    #if lenGap > 1:
    convertedScores = [getProbQuality(score[i]) for i in range(0,lenGap)]
    score = lambda sc : (92.21/100 * sc) + (4.79/100 * (1 - sc))
    #score = lambda sc : ((3.12/100 * sc) + (96.88/100 * (1 - sc)) + (4.79/100 * sc) + (95.21/100 * (1 - sc)))
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
@Definition: Calculation of overlap score for the overlap region
@Input parameters: Key (Seqs IDs), overlap data, scoring type
@Output parameters: overlap score
"""

def getOverlapScore(key, readData):
    #scoreList = [0]
    probabilityOverall = 0
    probabilityGaps = 0
    probabilityMatches = 0
    prob = 0
    L = 0
    seq1 = readData[0]
    seq2 = readData[5]
    score1 = readData[1]
    score2 = readData[6]
    r1S = readData[3]
    r1E = readData[4]
    r2S = readData[8]
    r2E = readData[9]
    result = []
    #alpha1 = (readData[4] - readData[3])/(readData[2] - readData[3] + readData[8])
    #alpha2 = (readData[9] - readData[8])/(readData[9] + readData[2] - readData[4])
    #alpha = max(alpha1, alpha2)
    k = max((readData[4] - readData[3]),(readData[9] - readData[8]))
    newAlpha = k/(k + min(readData[3],readData[8]) + min((readData[2] - readData[4]),(readData[7] - readData[9])))
    #c = 0
    #cig = readData[8].count("M") + readData[8].count("I") + readData[8].count("D")
    #print("CIGARDetails:",cig, readData[8].count("M"), readData[8].count("I"), readData[8].count("D"))
    #print("CIGAR:",readData[8])
    #print(len(seq1),len(seq2))
    stopPos = min((r1E-r1S),(r2E-r2S))
    r1 = r1S
    r2 = r2S
    cntr = 0
    print(stopPos,r1,r2,r1S,r1E,r2S,r2E)
    while L < stopPos:
        #print(cntr,r1,r2)
        probabilityBase = 0
        r1 = r1 + 1
        r2 = r2 + 1
        #print(r1,r2)
        for n in nt:
            try:
                sc = (probabilityQ(n,seq1[r1],score1[r1]) * probabilityQ(n,seq2[r2],score2[r2]))
                probabilityBase = probabilityBase + sc
            except:
                continue
        prob = prob + np.log(probabilityBase)
        L = L + 1
        #cntr = cntr + 1

    try:
        overlapScore = (np.exp(prob) ** (1/L)) * newAlpha
        result.append(overlapScore)
    except:
        result.append("ERROR")
    #overlapScore = (np.exp(prob) ** (1/L))
    result.append(overlapScore)
    #print(result)
    #result[3:3] = [alpha1,alpha2]
    result[1:1] = [newAlpha]
    return(result)
"""
    for j in range(0,3):
        #print(j)
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
                if j == 1:
                    pos1 = pos1 + num
                else:
                    sc = getGapRegionScore(tempScore1, num, scoreType)
                    prob = prob + np.log(sc)
                    pos1 = pos1 + num
                    L = L + 1
            elif char == "D":
                if j == 1:
                    pos2 = pos2 + num
                else:
                    sc = getGapRegionScore(tempScore2, num, scoreType)
                    prob = prob + np.log(sc)
                    pos2 = pos2 + num
                    L = L + 1
            elif char == "M":
                if j == 0:
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
        #print(alpha1,alpha2,alpha)

        try:
            overlapScore = (np.exp(prob) ** (1/L)) * newAlpha
            result.append(overlapScore)
        except:
            result.append("ERROR")
            continue
        #overlapScore = (np.exp(prob) ** (1/L))
        result.append(overlapScore)
        #print(result)
    #result[3:3] = [alpha1,alpha2]
    result[3:3] = [newAlpha]
    return(result)
"""

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
         fastqTemp[fastqData[i].split(" ")[0].strip("@").rstrip("\n")] = [fastqData[i+1].rstrip("\n"),fastqData[i+3].rstrip("\n")]
    i = i + 1

statFileName = args.output + "_stats_NoGaps.csv"
statFile = open(statFileName,"w+")
statFile.write("KEY \t READ1_LENGTH \t READ1_START \t READ1_END \t READ1_OVL_LEN \t READ2_LENGTH \t READ2_START \t READ2_END \t READ2_OVL_LEN \t MATCHES \t INSERTIONS \t DELETIONS \t GAPS \t MISMATCHES \n")
for overlapPair in pafData:
    tempData = list()
    ovl = overlapPair.split("\t")

    """
    if ovl[0] + "-" + ovl[5] == "origSeq_189-origSeq_190":

        #print(ovl)
        #print(fastqTemp[ovl[0]][0])
        #print(fastqTemp[ovl[0]][1])
        #print(fastqTemp[ovl[5]][0])
        #print(fastqTemp[ovl[5]][1])

        print(ovl)
        #print("Init:",ovl[20].split(":")[2].strip("\n"))
    """
    cigarIndex = [ovl.index(i) for i in ovl if i.startswith("cg")]
    cig = getExpandedCigar(ovl[cigarIndex[0]].split(":")[2].strip("\n"))
    gaps = cig.count("I") + cig.count("D")
    subsIndex = [ovl.index(i) for i in ovl if i.startswith("NM")]
    subsNumber = int(ovl[subsIndex[0]].split(":")[2]) - gaps
    statFile.write(str(ovl[0] + "-" + ovl[5]) + "\t" + str(ovl[1]) + "\t" + str(ovl[2]) + "\t" + str(ovl[3]) + "\t" + str(int(ovl[3]) - int(ovl[2])) + "\t" + str(ovl[6]) + "\t" + str(ovl[7]) + "\t" + str(ovl[8]) + "\t" + str(int(ovl[8]) - int(ovl[7])) + "\t" + str(cig.count("M")) + "\t" + str(cig.count("I")) + "\t" + str(cig.count("D")) + "\t" + str(gaps) + "\t" + str(subsNumber) + "\n")
    readPairData[ovl[0] + ":" + ovl[5]] = [fastqTemp[ovl[0]][0],fastqTemp[ovl[0]][1], int(ovl[1]), int(ovl[2]), int(ovl[3]), fastqTemp[ovl[5]][0], fastqTemp[ovl[5]][1], int(ovl[6]), int(ovl[7]),int(ovl[8]),ovl[cigarIndex[0]].split(":")[2].strip("\n")]
    #tempData = [list(fastq[ovl[0]].seq), fastq[ovl[0]].letter_annotations["phred_quality"], int(ovl[2]), int(ovl[3]), list(fastq[ovl[5]].seq),fastq[ovl[5]].letter_annotations["phred_quality"],int(ovl[7]), int(ovl[8]),ovl[20].split(":")[2].strip()]
    #readPairData[ovl[0] + "-" + ovl[5]] = [fastqTemp[ovl]]
statFile.close()

print(len(readPairData))
#c = 0
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

    scoreFileName = args.output + name + "_scores_NoGaps.csv"
    outputFile = open(scoreFileName,"w+")
    #outputFile.write("KEY \t ALPHA1 \t ALPHA2 \t GAPS \t MATCHES \t ALL \t OVERLAP_TYPE \n")
    outputFile.write("KEY \t ALPHA \t ALL \t OVERLAP_TYPE \n")

    for key in readPairData:
        #if key == "origSeq_125-origSeq_147":

            keyElements = key.split(":")
            if keyElements[0].split("_")[0] == keyElements[1].split("_")[0]:
                ovlType = "Good Overlap"
            else:
                ovlType = "Bad Overlap"

            results = getOverlapScore(key, readPairData[key])

            print(key,str(results[0]))

            #outputFile.write(key + "\t" + str(results[3]) + "\t" + str(results[4]) + "\t" + str(results[0]) + "\t" + str(results[1]) + "\t" + str(results[2]) + "\t" + ovlType + "\n")
            outputFile.write(key + "\t" + str(results[1]) + "\t" + str(results[0]) + "\t" + ovlType + "\n")
    outputFile.close()
