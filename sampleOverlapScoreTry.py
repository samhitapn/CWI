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
import re



nt = ["A","T","C","G"] # Bases
ntPattern = re.compile("[ATGC]")
gapPattern = re.compile("[-]")

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


"""
@Definition:
@Input parameters:
@Output parameters: Probability
"""
def findOverlapRegion (read1, read2):
    read1_start = min(read1.find("A"),read1.find("T"),read1.find("G"),read1.find("C"))
    read2_start = min(read2.find("A"),read2.find("T"),read2.find("G"),read2.find("C"))
    read1_end = max(read1.rfind("A"),read1.rfind("T"),read1.rfind("G"),read1.rfind("C"))
    read2_end = max(read2.rfind("A"),read2.rfind("T"),read2.rfind("G"),read2.rfind("C"))
    print(read1_start, read1_end, read2_start, read2_end)
    overlapStart = max(read1_start, read2_start)
    overlapEnd = min(read1_end, read2_end)
    return(overlapStart, overlapEnd)


"""
@Definition:
@Input parameters:
@Output parameters: Probability
"""
def reorderScores (read, score):
    i = 0
    while i < len(read):
        if read[i] == "-":
            score.insert(i,"-")
        i = i + 1
    return(score)

"""
@Definition:
@Input parameters:
@Output parameters: Probability
"""
def getGapRegionScore (read, score, pos):
    start = pos
    end = ntPattern.search(read,start).start()
    #print(start,end)
    lenGap = end - start
    #print(start, end, lenGap)
    if lenGap > 1:
        probSum = 0
        while start < end:
            #print(start,end,lenGap)
            probSum = probSum + (10/13 * getProbQuality(float(score[start]))) + (3/13 * (1 - getProbQuality(float(score[start]))))
            start = start + 1
        probGap = probSum/lenGap
    else:
        #print("ELSE")
        #print(start,end,lenGap)
        probGap = (10/13 * getProbQuality(float(score[start]))) + (3/13 * (1 - getProbQuality(float(score[start]))))
    return(probGap, end)



"""
@Definition:
@Input parameters:
@Output parameters:
"""
def overalpScoreCalculation(seqDetails):
    # Initiating required values
    probabilityOverall = 1

    # Getting the sequence and scores
    seqRead1 =  seqDetails.split(" ")[0]
    scoreRead1 = seqDetails.split(" ")[1].split(",")
    seqRead2 = seqDetails.split(" ")[2]
    scoreRead2 = seqDetails.split(" ")[3].split(",")

    # reordering the score positions based on alignment
    scoreRead1 = reorderScores(seqRead1,scoreRead1)
    scoreRead2 = reorderScores(seqRead2,scoreRead2)

    # Finding the overlap region
    overlapRegion = findOverlapRegion(seqRead1, seqRead2)
    startOverlap = overlapRegion[0]
    endOverlap = overlapRegion[1]
    L = endOverlap - startOverlap
    print(startOverlap, endOverlap, L)

    #outFile = open("scoreData.txt","w+")
    while startOverlap <= endOverlap:
        probabilityBase = 0

        #New code -> To include indels Option1
        """
        if seqRead1[startOverlap] == "-":
            print(startOverlap, seqRead1[startOverlap])
            end = ntPattern.search(seqRead1,startOverlap).start()
            print(end)
            startOverlap = end
            print(startOverlap)
            print("&&&&&&&")
        else:
            print(startOverlap, seqRead1[startOverlap])
            print("&&&&&&&")
            startOverlap = startOverlap + 1
        """
        # Gap in first read -> calculation based on read 2
        if seqRead1[startOverlap] == "-":
            pl = 1
            #print("BEFORE")
            #print(pl,startOverlap)
            #print(scoreRead2[startOverlap])
            #probabilityBase = (3/13 * getProbQuality(float(scoreRead2[startOverlap]))) + (10/13 * (1 - getProbQuality(float(scoreRead2[startOverlap]))))
            #probabilityBase = (10/13 * float(scoreRead2[startOverlap])) + (3/13 * (1 - float(scoreRead2[startOverlap])))
            gapDetails = getGapRegionScore(seqRead1,scoreRead2,startOverlap)
            probabilityBase = gapDetails[0]
            startOverlap = gapDetails[1]

        # Gap in second read -> calculation based on read 1
        elif seqRead2[startOverlap] == "-":
            pl = 2
            #print("BEFORE")
            #print(pl,startOverlap)
            #print(scoreRead1[startOverlap])
            #probabilityBase = (10/13 * float(scoreRead1[startOverlap])) + (3/13 * (1 - float(scoreRead1[startOverlap])))
            #end = ntPattern.search(seqRead2,startOverlap).start()
            #print(end)
            #gapDetails = getGapRegionScore(seqRead1,scoreRead1,startOverlap)
            gapDetails = getGapRegionScore(seqRead2,scoreRead1,startOverlap)
            probabilityBase = gapDetails[0]
            startOverlap = gapDetails[1]
            #print("AFTER")
            #print(pl, startOverlap)



        # Existing score calculation
        else:
            pl = 3
            #print("BEFORE")
            #print(pl,startOverlap)
            for n in nt:
               #print(n)
               #probabilityBase = probabilityBase + (probabilityQ(n,seqRead1[startOverlap],getProbQuality(float(scoreRead1[startOverlap]))) * probabilityQ(n,seqRead2[startOverlap],getProbQuality(float(scoreRead2[startOverlap]))))
               probabilityBase = probabilityBase + (probabilityQ(n,seqRead1[startOverlap],float(scoreRead1[startOverlap])) * probabilityQ(n,seqRead2[startOverlap],float(scoreRead2[startOverlap])))
            startOverlap = startOverlap + 1
            #print("AFTER")
            #print(pl, startOverlap)
        #print("$$$$$$$$$$$$$$$$$$$$$$$$$$")
            #print("I am in %f", pl)
        #print(probabilityBase)
        probabilityOverall = probabilityOverall * probabilityBase

        print(startOverlap-1, seqRead1[startOverlap-1], seqRead2[startOverlap-1], scoreRead1[startOverlap-1] , scoreRead2[startOverlap-1] , probabilityBase , probabilityOverall , pl)

      # Overlap score
    overlapScore = probabilityOverall ** 1/L
    #print(overlapScore)
    #outFile.close()
    return (overlapScore)

# MAIN
sequences = sio.to_dict(sio.parse("data/Sample_AllReads.fastq","fastq"))
print("^^^^^^SEQUENCES READ^^^^^^\n")

overlapCount = 1
alignment = {}
# Get the overlap pairs and details from the PAF files
with open("data/Sample_AllReads_Overlaps.paf","r") as f:
    for ovlPair in f.readlines():
        if overlapCount <= 1:
        #print(ovlPair)
        #print("^^^^^^^^^^^^^^^^^^^^")

            ovl = ovlPair.split("\t")
            if ovl[0] != ovl[5]:
                #print(overlapCount)
                start = time()
                al = pairwise2.align.globalxx(sequences[ovl[0]].seq,sequences[ovl[5]].seq)
                stop = time()
                print(stop - start)
                alignment[ovl[0] + ' ' + ovl[5]] = al[0][0] + " " + str(sequences[ovl[0]].letter_annotations["phred_quality"]).strip('[]').replace(" ","") + " " + al[0][1] + " " + str(sequences[ovl[0]].letter_annotations["phred_quality"]).strip('[]').replace(" ","")
                overlapCount = overlapCount + 1
                #print(al[0][0])
                #print(al[0][1])
                #print("$$$$$$$$$$$$$$$$$$")
        else:
            break
#print(alignment)

# Getting the scores for each overlap pair
for key in alignment:
        #print(alignment[key].split(" ")[0])
        #print(alignment[key].split(" ")[2])
        print(key)
        overlapScore = overalpScoreCalculation(alignment[key])
        alignment[key] = alignment[key] + " " + str(overlapScore)
        print(alignment[key].split(" ")[4])
        print("\n ############ \n")
