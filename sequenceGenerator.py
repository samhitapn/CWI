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
    seqList = re.findall("[ATGC]",originalSequence)
    print(i)
    mutationPercent = i/100 * 100000
    #seed = 123
    print(mutationPercent)
    pos = random.sample(range(0,100000),int(mutationPercent))
    for n in tqdm(pos):
        nReplace = [j for j in nt if j != originalSequence[n]]
        seqList[n] = random.choice(nReplace)
    seq = "".join(seqList)
    print(len(seq))
    #print(seq)
    fileName = "seq" + str(i)
    with open("data/" + fileName + ".fasta", "w") as file:
            file.write(">" + fileName + "\n" + seq)
    print("Written")
