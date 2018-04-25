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
originalSequence = rstr.rstr('ATGC', 10000000)
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

for i in tqdm(range(1,2)):
    seq = originalSequence
    print(i)
    mutationPercent = i/100 * 10000000
    seed = 123
    print(mutationPercent)
    pos = random.sample(range(0,10000000),int(mutationPercent))
    for n in tqdm(pos):
        nReplace = [j for j in nt if j != originalSequence[n]]
        seq = seq.replace(seq[n],random.choice(nReplace))
    #print(seq)
    fileName = "seq" + str(i)
    with open("data/" + fileName + ".fasta", "w") as file:
            file.write(">" + fileName + "\n" + seq)
    print("Written")
