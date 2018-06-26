"""
@created: 24.June.2018
@author: samhitapn
@purpose: Generate and analyse sequencing error-free reads
@modified: 26.June.2018
"""

import re
import os
import argparse

cigarPattern = re.compile('([0-9]*)([DMSH])')

"""
@Definition: Converting the CIGAR string to the expanded form
@Input parameters: CIGAR string
@Output parameters: expanded CIGAR string
"""
def getExpandedCigar(cigar):
    cigarSeq = []
    start = True
    for num, char in cigarPattern.findall(cigar):
        if num:
            num = int(num)
        else:
            num = 1
        if start and (char == "S" or char == "H"):
            s = num
            start = False
        cigarSeq.append("".join(char * num))
    cigarSeq = "".join(cigarSeq)
    #reqLength = cigarSeq.count("S") + cigarSeq.count("H") + cigarSeq.count("D") + cigarSeq.count("M")
    return(cigarSeq,s)


# Arguement parser
parser = argparse.ArgumentParser(description='Parser for input files and output files')
parser.add_argument('-f','--file', help='File Name',required=True)
args = parser.parse_args()

# Input data
fastqTemp = dict()

with open(args.file + ".sam") as s:
    sam = s.readlines()
with open(args.file + ".fastq") as fq:
    fastqData = fq.readlines()
for i in range(0,len(fastqData)):
    if fastqData[i].startswith("@"):
         fastqTemp[fastqData[i].split(" ")[0].strip("@")] = [fastqData[i+1].rstrip("\n"),fastqData[i+3].rstrip("\n")]
    i = i + 1

# Replacing fastq Dictionary
c = 0
for line in sam:
    if c >= 2:
        temp = line.split("\t")
        readName = temp[0]
        pos = getExpandedCigar(temp[5])
        reqLength = pos[0].count("S") + pos[0].count("H") + pos[0].count("D") + pos[0].count("M")
        startPos = temp[3] - (1 + reqLength)
        endPos = startPos + pos[0]
        with open("../fasta_100000_indel/" + args.file + ".fasta") as fasta:
            seq = fasta.readlines()
        sequence = seq[1][startPos:endPos]
        fastqTemp[readName][0] = sequence

# Writing replaced error-free fastq file
fileNew_fastq = args.file + "_errorFree.fastq"
with open(fileNew_fastq,"w+") as fw:
    for key in fastqTemp:
        fw.write("@" + key + "\n" + fastqTemp[key][0] + "\n+\n" + fastqTemp[key][1] + "\n")

# Generating overlaps from error-free reads
for i in ["eb0","eb10","eb100","eb1000"]:
    fileNew_paf = i + "/" + args.file + "_errorFree.paf"
    cmd = "minimap2 -x ava-ont " + fileNew_fastq + " " + fileNew_fastq + " -c --end-bonus " + i[2:] + " > " + fileNew_paf
    os.system(cmd)

# Parsing CIGAR string from both PAF files for gaps
for i in ["eb0","eb10","eb100","eb1000"]:
    print(i)
    os.system("pwd")
    os.chdir(i)
    file_paf = args.file + ".paf"
    with open(file_paf) as paf:
        pafData = paf.readlines()
    with open(fileNew_paf) as pafNew:
        pafData_New = pafNew.readline()
    print("PAF DATA RECEIVED")
    with open(args.file + "_CIGAR.csv","w+") as oldCigar:
        oldCigar.write("KEY \t GAPS \t MATCHES \t DELETIONS \t INSERTIONS")
        for pair in pafData:
            ovl = pair.split("\t")
            cigarIndex = [ovl.index(i) for i in ovl if i.startswith("cg")]
            expCigar = getExpandedCigar(ovl[cigarIndex[0]].split(":")[2].strip("\n"))
            gaps = expCigar[0].count("D") + expCigar[0].count("I")
            oldCigar.write(str(ovl[0] + "-" + ovl[5]) + "\t" + str(gaps) + "\t" + str(expCigar[0].count("M")) + "\t" + expCigar[0].count("D") + "\t" + expCigar[0].count("I"))

    with open(args.file + "_errorFree_CIGAR.csv","w+") as newCigar:
        newCigar.write("KEY \t GAPS \t MATCHES \t DELETIONS \t INSERTIONS")
        for pair in pafData_New:
            ovl = pair.split("\t")
            cigarIndex = [ovl.index(i) for i in ovl if i.startswith("cg")]
            expCigar = getExpandedCigar(ovl[cigarIndex[0]].split(":")[2].strip("\n"))
            gaps = expCigar[0].count("D") + expCigar[0].count("I")
            newCigar.write(str(ovl[0] + "-" + ovl[5]) + "\t" + str(gaps) + "\t" + str(expCigar[0].count("M")) + "\t" + expCigar[0].count("D") + "\t" + expCigar[0].count("I"))
    print("PAF STATS SAVED")
    os.chdir("..")
