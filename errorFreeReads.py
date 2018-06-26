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
    clipCount = 0
    for num, char in cigarPattern.findall(cigar):
        if num:
            num = int(num)
        else:
            num = 1
        if start and (char == "S" or char == "H"):
            clipCount = num
            start = False
        cigarSeq.append("".join(char * num))
    cigarSeq = "".join(cigarSeq)
    #reqLength = cigarSeq.count("S") + cigarSeq.count("H") + cigarSeq.count("D") + cigarSeq.count("M")
    return(cigarSeq,clipCount)


# Arguement parser
"""
parser = argparse.ArgumentParser(description='Parser for input files and output files')
parser.add_argument('-f','--file', help='File Name',required=True)
args = parser.parse_args()
"""
# Input data
for file in ["origSeq", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6", "seq7", "seq8", "seq9", "seq10"]:
    fastqTemp = dict()

    with open(file + ".sam") as s:
        sam = s.readlines()
    with open(file + ".fastq") as fq:
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
            with open("../fasta_100000_indel/" + file + ".fasta") as fasta:
                seq = fasta.readlines()
            sequence = seq[1][startPos:endPos]
            fastqTemp[readName][0] = sequence

    # Writing replaced error-free fastq file
    fileNew_fastq = file + "_errorFree.fastq"
    with open(fileNew_fastq,"w+") as fw:
        for key in fastqTemp:
            fw.write("@" + key + "\n" + fastqTemp[key][0] + "\n+\n" + fastqTemp[key][1] + "\n")

# Concatenating all files

cmdCat = "for f in *_errorFree.fastq; do (cat \"${f}\"; echo) >> all/allMerged_200_errorFree.fastq; done"
os.system(cmdCat)

# Generating overlaps from error-free reads
os.chdir("all")
for j in ["EB0","EB10","EB100","EB1000"]:
    file_paf = "allMerged_" + j + "_errorFree.paf"
    cmd = "minimap2 -x ava-ont " + fileNew_fastq + " " + fileNew_fastq + " -c --end-bonus " + j[2:] + " > " + file_paf
    os.system(cmd)

# Parsing CIGAR string from both PAF files for gaps
#fileList = [f for f in os.listdir(".") if "allMerged_errorFree_" in f]
for i in ["EB0","EB10","EB100","EB1000"]:
    print(i)
    #os.system("pwd")
    #os.chdir(i)
    file_paf = "allMerged_" + i + ".paf"
    with open(file_paf) as paf:
        pafData = paf.readlines()
    fileNew_paf = "allMerged_" + i + "_errorFree.paf"
    with open(fileNew_paf) as pafNew:
        pafData_New = pafNew.readline()
    print("PAF DATA RECEIVED")
    with open(i + "_CIGAR.csv","w+") as oldCigar:
        oldCigar.write("KEY \t GAPS \t MATCHES \t DELETIONS \t INSERTIONS")
        for oldPair in pafData:
            ovl = oldPair.split("\t")
            cigarIndex = [ovl.index(of) for of in ovl if of.startswith("cg")]
            expCigar = getExpandedCigar(ovl[cigarIndex[0]].split(":")[2].strip("\n"))
            gaps = expCigar[0].count("D") + expCigar[0].count("I")
            oldCigar.write(str(ovl[0] + "-" + ovl[5]) + "\t" + str(gaps) + "\t" + str(expCigar[0].count("M")) + "\t" + str(expCigar[0].count("D")) + "\t" + str(expCigar[0].count("I")))
    print("OLD DONE")
    with open(i + "_errorFree_CIGAR.csv","w+") as newCigar:
        newCigar.write("KEY \t GAPS \t MATCHES \t DELETIONS \t INSERTIONS")
        for newPair in pafData_New:
            ovl = newPair.split("\t")
            cigarIndex = [ovl.index(nf) for nf in ovl if nf.startswith("cg")]
            print(cigarIndex)
            expCigar = getExpandedCigar(ovl[cigarIndex[0]].split(":")[2].strip("\n"))
            gaps = expCigar[0].count("D") + expCigar[0].count("I")
            newCigar.write(str(ovl[0] + "-" + ovl[5]) + "\t" + str(gaps) + "\t" + str(expCigar[0].count("M")) + "\t" + str(expCigar[0].count("D")) + "\t" + str(expCigar[0].count("I")))
    print("NEW DONE")
    print("PAF STATS SAVED")
    os.chdir("..")
