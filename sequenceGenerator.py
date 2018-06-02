#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created: 24.April.2018
@author: samhitapn
@purpose: Random Sequence generator
@modified: 24.April.2018
"""

from Bio import pairwise2
from Bio import SeqIO as sio
import numpy as np
import math as mt
import itertools
from time import time
import re
import random
import rstr
from tqdm import tqdm

nt = ["A","T","G","C"]

startGenerationOrig = time()
originalSequence = rstr.rstr('ATGC', 100000)
stopGenerationOrig = time()
print(stopGenerationOrig - startGenerationOrig)
print("Reference Done")
startOrigWrite = time()
with open("data/origSeq.fasta","w") as origFile:
    origFile.write(">OriginalSequenceReference\n" + originalSequence)
    #origFile.write(originalSequence)
stopOrigWrite = time()
print(stopOrigWrite - startOrigWrite)
print("Written")

for i in tqdm(range(1,11)):
    print(i)
    seqList = re.findall("[ATGC]",originalSequence)
    print("INITIAL",len(seqList))
    #print(i)
    mutationNumber = i/100 * 100000

    # DELETIONS
    delNumber = 3.2/100 * mutationNumber
    #pos = random.sample(range(0,100000),int(mutationNumber))
    # CHECKING THE VALIDITY OF THE POSITIONS
    while repeat:
        delPos = random.sample(range(0,len(seqList)),int(delNumber))
        delPos = sort(delPos)
        delChoices = random.choices([1,2,3,4,5],k = int(delNumber))
        repeat = FALSE
        for j in range(0,delNumber - 1):
             temp = delPos[j] + delChoices[j]
             if temp > delPos[j + 1]:
                 repeat = TRUE
                 break
    # REPLACING THE BASES TO BE DELETED
    for d in range(0,delNumber):
        start = delPos[d]
        end = delPos[d] + choices[d]
        seqList[delPos[d]:delPos[d] + choices[d]] = "X"
    # DELETING THE BASES
    seqList = [b for b in seqList if b != "X"]
    print("AFTER DELETION",len(seqList))


    # SUBSTITUIONS
    subNumber = 93.2/100 * mutationNumber
    #pos2 = [x for x in pos1 if x not in insPos]
    subPos = random.sample(range(0,len(seqList)), int(subNumber))
    for n in subPos:
        nReplace = [b for b in nt if b != originalSequence[n]]
        seqList[n] = random.choice(nReplace)
    print("AFTER SUBSTITUION",len(seqList))

    #INSERTIONS
    insNumber = 3.6/100 * mutationNumber
    #pos1 = [x for x in pos if x not in delPos]
    insPos = random.sample(range(0,len(seqList) - 1),int(insNumber))
    # CHECKING THE VALIDITY OF THE POSITIONS
    while repeat:
        insPos = random.sample(range(0,len(seqList)),int(insNumber))
        insPos = sort(insPos)
        insChoices = random.choices([1,2,3,4,5],k = int(insNumber))
        repeat = FALSE
        for k in range(0,insNumber - 1):
             temp = inSPos[k] + insChoices[k]
             if temp > insPos[k + 1]:
                 repeat = TRUE
                 break
    # INSERTING THE REQUIRED BASES
    counter = 0
    for d in range(0,insNumber):
        insertions = random.choices(["A","T","G","C"],k = insChoices[d])
        seqList[insPos[d] + counter : insPos[d] + counter] = insertions
        counter = counter + insChoices[d]
    print("AFTER INSERTION",len(seqList))

    """
    # JOIN THE SEQUENCE LIST
    seq = "".join(seqList)
    print(len(seq))
    #print(seq)
    fileName = "seq" + str(i)
    with open("data/" + fileName + ".fasta", "w") as file:
            file.write(">" + fileName + "\n" + seq)
    print("Written")
    """
