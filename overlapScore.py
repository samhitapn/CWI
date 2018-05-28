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
def getGapRegionScore (score, lenGap):
    ovlScore = 0
    #print(scoreList)
    #scoreList = []
    #if lenGap > 1:

    for i in range(0, lenGap):
        #print(score[i])
        #score[i] = getProbQuality(score[i])
        #scoreList = scoreList.append(getProbQuality(score[i]))
        #print("Score&&&&&",ord(score[i]))
        ovlScore = ovlScore + getProbQuality(score[i])
    ovlScore = ovlScore/lenGap
    prob = (10/13 * ovlScore) + (3/13 * (1 - ovlScore))
    #print("PROB", prob)
    #print(scoreList)
    return(prob)

"""
@Definition: Mapping the sequences, getting the sequence and the scores and converting the CIGAR string to the alignment
@Input parameters:
@Output parameters:
"""

def getOverlapScore(key, readData):
    #scoreList = [0]
    probabilityOverall = 0
    probabilityGaps = 0
    probabilityMatches = 0
    prob = 0
    seq1 = readData[0]
    seq2 = readData[4]
    score1 = readData[1]
    score2 = readData[5]
    pos1 = readData[2]
    pos2 = readData[6]
    L = 0
    errorDet = list()
    #c = 0
    #cig = readData[8].count("M") + readData[8].count("I") + readData[8].count("D")
    #print("CIGARDetails:",cig, readData[8].count("M"), readData[8].count("I"), readData[8].count("D"))
    #print("CIGAR:",readData[8])
    #print(len(seq1),len(seq2))
    for num, char in cigarPattern.findall(readData[8]):
        try:
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
                #print("!!!!",char, num, "*******")
                sc = getGapRegionScore(tempScore1, num)
                prob = prob + np.log(sc)
                #print(sc[1])
                #scoreList = sc[1]
                #assert 0 <= probabilityOverall <= 1, print(char, pos1, pos2, probabilityOverall)
                #assert probabilityOverall <= 0.000001, print(num,char, pos1, pos2, probabilityOverall,tempScore2,tempSeq2)
                pos1 = pos1 + num
                L = L + 1
                #print("I:",num,pos1,pos2)
                #print("@@@@@ I-range",num, probabilityOverall, tempScore1,tempScore2, pos1, pos2, L, c)

            elif char == "D":
                #print("!!!!",char, num, "*******")
                sc = getGapRegionScore(tempScore2, num)
                prob = prob + np.log(sc)
                #print(sc[1])
                #scoreList = sc[1]
                #probabilityOverall = probabilityOverall * getGapRegionScore(tempScore1, num)
                #assert 0 <= probabilityOverall <= 1, print(char, pos1, pos2, probabilityOverall)
                #assert probabilityOverall <= 0.000001, print(num,char, pos1, pos2, probabilityOverall,tempScore1,tempSeq1)
                pos2 = pos2 + num
                L = L + 1
                #print("D:",num,pos1,pos2)
                    #print("@@@@@ D-range" , num,probabilityOverall, tempScore1,tempScore2, pos1, pos2, L, c)

            elif char == "M":

            #if char == "M":
                #print("!!!!", char, num, "*******")

                for i in range(0, num):
                    probabilityBase = 0
                    #print(tempScore1[i],getProbQuality(tempScore1[i]),tempScore2[i],getProbQuality(tempScore2[i]))
                    for n in nt:
                        sc = (probabilityQ(n,tempSeq1[i],tempScore1[i]) * probabilityQ(n,tempSeq2[i],tempScore2[i]))
                        probabilityBase = probabilityBase + sc
                        #assert 0 <= probabilityBase <= 1, print(char, pos1, pos2, probabilityBase, "Base")
                        #print(n, sc, probabilityBase, tempSeq1[i],ord(tempScore1[i]),tempSeq2[i],ord(tempScore2[i]))
                        #probabilityBase = probabilityBase + (probabilityQ(n,tempSeq1[i],np.float128(tempScore1[i])) * probabilityQ(n,tempSeq2[i],np.float128(tempScore2[i])))
                    prob = prob + np.log(probabilityBase)
                    #assert 0 <= probabilityOverall <= 1, print(char, pos1, pos2, probabilityOverall)
                    L = L + 1
                    #print("@@@@@ M-location", probabilityOverall,L,c,pos1,pos2)
                #print("@@@@@ M-range", num,probabilityOverall, tempScore1,tempScore2, pos1,pos2,L,c)

                pos1 = pos1 + num
                pos2 = pos2 + num
                #print("M:",num,pos1,pos2)

        except (IndexError, UnboundLocalError):
            #print("ERROR",probabilityOverall)
            #print("TTTTTT",probabilityOverall,np.exp(probabilityOverall))
            overlapScore = np.exp(prob) ** (1/L)
            result = [1, overlapScore]
            continue

        else:
            #print("TTTTTT",probabilityOverall,np.exp(probabilityOverall))
            overlapScore = np.exp(prob) ** (1/L)
            #print(L, probabilityOverall,np.exp(probabilityOverall),overlapScore)
            result = [0, overlapScore]
        #brEak

    return(result)


# MAIN
# Arguement parsing
parser = argparse.ArgumentParser(description='Parser for input files and output files')
parser.add_argument('-p','--paf', help='PAF file name',required=True)
parser.add_argument('-f','--fastq',help='Fastq file name', required=True)
parser.add_argument('-o','--output',help='Fastq file name', required=True)
#parser.add_argument('-g','--gap',help='Gap region method max (0), min (1), avg (2), all(3)', required=True, type = int)
#parser.add_argument('-t','--scoreArea',help='Score only gaps and/or substitutions: Gaps(0), Substitutions(1), Both(2), All(3)', required=True, type = int)
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
statFile.write("KEY \t READ1_START \t READ1_END \t READ1_OVL_LEN \t READ2_START \t READ2_END \t READ2_OVL_LEN \t MATCHES \t INSERTIONS \t DELETIONS \n")
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
    statFile.write(str(ovl[0] + "-" + ovl[5]) + "\t" + str(ovl[2]) + "\t" + str(ovl[3]) + "\t" + str(int(ovl[3]) - int(ovl[2])) + "\t" + str(ovl[7]) + "\t" + str(ovl[8]) + "\t" + str(int(ovl[8]) - int(ovl[7])) + "\t" + str(cig.count("M")) + "\t" + str(cig.count("I")) + "\t" + str(cig.count("D")) + "\n")
    readPairData[ovl[0] + "-" + ovl[5]] = [fastqTemp[ovl[0]][0],fastqTemp[ovl[0]][1],int(ovl[2]), int(ovl[3]),fastqTemp[ovl[5]][0],fastqTemp[ovl[5]][1],int(ovl[7]),int(ovl[8]),ovl[cigarIndex[0]].split(":")[2].strip("\n")]
    #tempData = [list(fastq[ovl[0]].seq), fastq[ovl[0]].letter_annotations["phred_quality"], int(ovl[2]), int(ovl[3]), list(fastq[ovl[5]].seq),fastq[ovl[5]].letter_annotations["phred_quality"],int(ovl[7]), int(ovl[8]),ovl[20].split(":")[2].strip()]
    #readPairData[ovl[0] + "-" + ovl[5]] = [fastqTemp[ovl]]
statFile.close()

print(len(readPairData))
c = 0
scoreFileName = args.output + "_scores.csv"
outputFile = open(scoreFileName,"w+")
for key in readPairData:
    c = c + 1
    #if c == 13931:
    print(key)
    #args.gap = [args.gap]
    #print(args.gap,type(args.gap))

    results = getOverlapScore(key, readPairData[key])
    keyElements = key.split("-")

    if keyElements[0].split("_")[0] == keyElements[1].split("_")[0]:
        ovlType = "Good Overlap"
    else:
        ovlType = "Bad Overlap"

    if results[2] == 1:
        error = "IndexError"
    else:
        error = "No error"

    print(c,key,str(results[3]))
    outputFile.write(key + "\t" + error + "\t" + str(results[3]) + "\t" + ovlType + "\n")
outputFile.close()
