"""
@created: 24.June.2018
@author: samhitapn
@purpose: Generate and analyse sequencing error-free reads
@modified: 26.June.2018
"""

import re
import os
import argparse

cigarPattern = re.compile('([0-9]*)([IDMSH])')

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

parser = argparse.ArgumentParser(description='Parser for input files and output files')
parser.add_argument('-f','--file', help='File Name',required=True)
parser.add_argument('-eb','--endbonus', help='end bonus choice',required=True)
#parser.add_argument('-o','--old ', help='File Name',required=True)
args = parser.parse_args()

files = args.file
files = files.split(",")
eb = args.endbonus
eb = eb.split(",")
# Input data

for file in files:
    print(file)
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
    #c = 0
    with open(file + ".fasta") as fasta:
        seq = fasta.readlines()
    c = 0
    for line in sam:
        print(c)
        if not line.startswith("@"):
            temp = line.split("\t")
            readName = temp[0]
            pos = getExpandedCigar(temp[5])
            reqLength = pos[0].count("S") + pos[0].count("H") + pos[0].count("D") + pos[0].count("M")
            trimLength = max(0, pos[1] - (int(temp[3]) - 1))
            startPos = max(0, int(temp[3]) - (1 + pos[1]))
            endPos = min(len(seq[1]), startPos + reqLength - trimLength)
            if startPos >= endPos:
                print(temp[:6])
                print(trimLength)
                print(reqLength)
                print(startPos, endPos)
            try:
                assert startPos < endPos
            except AssertionError:
                continue
            sequence = seq[1][startPos:endPos]
            if len(sequence) == 0:
                print(startPos, endPos)
            assert len(sequence) > 0
            fastqTemp[readName][0] = sequence
            #if readName == 'origSeq_192':
                #print(fastqTemp[readName][0])
            c = c + 1
    # Writing replaced error-free fastq file
    fileNew_fastq = file + "_errorFree.fastq"
    with open(fileNew_fastq,"w+") as fw:
        fw.seek(0)
        #n = "\n"
        for key in fastqTemp:
            #line = "@",key,n,fastqTemp[key][0],n,"+",n,fastqTemp[key][1],n
            #if key == 'origSeq_192':
                #print(fastqTemp[key][0])
            fw.writelines("@" + key + fastqTemp[key][0] + "\n+\n" + fastqTemp[key][1] + "\n")
            #fw.writelines(">" + key + "\n" + fastqTemp[key][0] + "\n")

    fileNew_fasta = file + "_errorFree.fasta"
    with open(fileNew_fasta,"w+") as faw:
        faw.seek(0)
        #n = "\n"
        for key in fastqTemp:
            #line = "@",key,n,fastqTemp[key][0],n,"+",n,fastqTemp[key][1],n
            #if key == 'origSeq_192':
                #print(fastqTemp[key][0])
            faw.writelines("@" + key + fastqTemp[key][0] + "\n")

print("DONE REPLACING INDIVIDUAL FASTQ  AND FASTA")

# Concatenating all files
#cmdCat = "for f in *_errorFree.fastq; do (cat \"${f}\"; echo) >> all/allMerged_200_errorFree.fastq; done"
mergeFile = ""
for file in files:
    mergeFile = mergeFile + file
cmdCat = "for f in *_errorFree.fasta; do (cat \"${f}\"; echo) >>" + mergeFile + "_errorFree.fasta; done"
os.system(cmdCat)
print("DONE COMBINING FASTQ's")


"""
# Generating overlaps from error-free reads
#os.chdir("all")
#os.system("pwd")
for i in eb:
    print(i)
    file_paf = mergeFile + i + "_error.paf"
    fileNew_paf = mergeFile + i + "_errorFree.paf"
    #cmdOld = "minimap2 -x ava-ont" +  mergeFile + ".fasta" +  mergeFile + ".fasta -c --end-bonus " + i[2:] + " > " + file_paf
    cmdNew = "minimap2 -x ava-ont " +  mergeFile + "_errorFree.fasta " +  mergeFile + "_errorFree.fasta -c --end-bonus " + i[2:] + " > " + fileNew_paf
    #print(cmd)
        #os.system(cmdOld)
    os.system(cmdNew)
    print("OVERLAP GENRATION DONE")
    #os.system("du -shx *|sort -rh")
    #with open(file_paf) as paf:
    #    pafData = paf.readlines()
    with open(fileNew_paf) as pafNew:
        pafData_New = pafNew.readlines()
    print("PAF DATA RECEIVED")
    with open(mergeFile + i + "_errorFree_CIGAR.csv","w+") as newCigar:
        newCigar.seek(0)
        newCigar.write("KEY \t GAPS \t MATCHES \t DELETIONS \t INSERTIONS \t MISMATCHES \n")
        for newPair in pafData_New:
            ovl = newPair.split("\t")
            cigarIndex = [ovl.index(nf) for nf in ovl if nf.startswith("cg")]
            expCigar = getExpandedCigar(ovl[cigarIndex[0]].split(":")[2].strip("\n"))
            gaps = expCigar[0].count("D") + expCigar[0].count("I")
            subsIndex = [ovl.index(nf) for nf in ovl if nf.startswith("NM")]
            subsNumber = int(ovl[subsIndex[0]].split(":")[2]) - gaps
            newCigar.write(str(ovl[0] + ":" + ovl[5]) + "\t" + str(gaps) + "\t" + str(expCigar[0].count("M")) + "\t" + str(expCigar[0].count("D")) + "\t" + str(expCigar[0].count("I")) + "\t" + str(subsNumber) + "\n")
    print("NEW DONE")
    print("PAF STATS SAVED")

#####################



# Parsing CIGAR string from both PAF files for gaps
for i in ["EB0","EB10","EB100","EB1000"]:
    print(i)
    #os.system("pwd")
    #os.chdir(i)
    file_paf = "allMerged_" + i + ".paf"
    with open(file_paf) as paf:
        pafData = paf.readlines()
    fileNew_paf = "allMerged_" + i + "_errorFree.paf"
    with open(fileNew_paf) as pafNew:
        pafData_New = pafNew.readlines()
    print("PAF DATA RECEIVED")
    with open(i + "_CIGAR.csv","w+") as oldCigar:
        oldCigar.seek(0)
        oldCigar.write("KEY \t GAPS \t MATCHES \t DELETIONS \t INSERTIONS \n")
        for oldPair in pafData:
            ovl = oldPair.split("\t")
            cigarIndex = [ovl.index(of) for of in ovl if of.startswith("cg")]
            expCigar = getExpandedCigar(ovl[cigarIndex[0]].split(":")[2].strip("\n"))
            gaps = expCigar[0].count("D") + expCigar[0].count("I")
            oldCigar.write(str(ovl[0] + "-" + ovl[5]) + "\t" + str(gaps) + "\t" + str(expCigar[0].count("M")) + "\t" + str(expCigar[0].count("D")) + "\t" + str(expCigar[0].count("I")) + "\n")
    print("OLD DONE")
    with open(i + "_errorFree_CIGAR.csv","w+") as newCigar:
        newCigar.seek(0)
        newCigar.write("KEY \t GAPS \t MATCHES \t DELETIONS \t INSERTIONS \n")
        for newPair in pafData_New:
            ovl = newPair.split("\t")
            cigarIndex = [ovl.index(nf) for nf in ovl if nf.startswith("cg")]
            expCigar = getExpandedCigar(ovl[cigarIndex[0]].split(":")[2].strip("\n"))
            gaps = expCigar[0].count("D") + expCigar[0].count("I")
            newCigar.write(str(ovl[0] + "-" + ovl[5]) + "\t" + str(gaps) + "\t" + str(expCigar[0].count("M")) + "\t" + str(expCigar[0].count("D")) + "\t" + str(expCigar[0].count("I")) + "\n")
    print("NEW DONE")
    print("PAF STATS SAVED")


        with open(i + "_CIGAR.csv","w+") as oldCigar:
            oldCigar.seek(0)
            oldCigar.write("KEY \t GAPS \t MATCHES \t DELETIONS \t INSERTIONS \t MISMATCHES \n")
            for oldPair in pafData:
                ovl = oldPair.split("\t")
                cigarIndex = [ovl.index(of) for of in ovl if of.startswith("cg")]
                expCigar = getExpandedCigar(ovl[cigarIndex[0]].split(":")[2].strip("\n"))
                gaps = expCigar[0].count("D") + expCigar[0].count("I")
                subsIndex = [ovl.index(of) for of in ovl if of.startswith("NM")]
                subsNumber = int(ovl[subsIndex[0]].split(":")[2]) - gaps
                oldCigar.write(str(ovl[0] + "-" + ovl[5]) + "\t" + str(gaps) + "\t" + str(expCigar[0].count("M")) + "\t" + str(expCigar[0].count("D")) + "\t" + str(expCigar[0].count("I")) + "\t" + str(subsNumber) + "\n")
        print("OLD DONE")
os.chdir("..")
"""
