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

nt = ["A","T","G","C"]

originalSequence = rstr.rstr('ATGC', 10000000)
with open("data/origSeq.fasta","w") as origFile:
    origFile.write(">OriginalSequenceReference\n" + originalSequence)
    #origFile.write(originalSequence)

for i in range(1,11):
    mutationPercent = i/100 * 10000000
    pos = random.sample(range(0,999999),int(mutationPercent))
    for n in pos:
        nReplace = [j for j in nt if j != originalSequence[n]]
        seq = originalSequence.replace(originalSequence[n],choice(nReplace))
        fileName = "seq" + str(i)
        with open("data/" + fileName + ".fasta") as file:
            file.write(">" + fileName + "\n" + seq)
