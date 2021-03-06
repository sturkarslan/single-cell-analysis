from __future__ import division
from subprocess import Popen, PIPE, STDOUT
import glob, sys, os, string, subprocess

folders = glob.glob("results-03/dvh/DvH_03_16-36162127/")
variantFile = "dvh-UA3-152-03-singlecell-variants-2callers-80percent-2cells_noan-bed.txt"
fastaFile = "reference/Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.fasta"

def runBamReadCount():
    program = "~/github/bam-readcount/bin/bam-readcount"
    parameters = "-q 0 -w 5"
    # # check the header of bam file and inspect which chromosome name is defined
    # cmd0 = "samtools idxstats %s | cut -f 1 | head -1" %(bamFile)
    # pipe = subprocess.Popen(cmd0, shell=True, stdout=PIPE).stdout
    # bamheader = pipe.read()
    # # assign different fasta file for different chromosome names
    # if bamheader == "NC_005863\n":
    #     fastaFile = "../dvh-coculture-rnaseq/dvh-single-cells/reference/Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.newname.fasta"
    #     # file to get list of all variants in all single cells
    #     variantFile = "dvh-UA3-152-03-singlecell-variants-2callers-80percent-2cells_noan-bed.txt"
    # else:
    #     fastaFile = "../dvh-coculture-rnaseq/dvh-single-cells/reference/Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.fasta"
    #     # file to get list of all variants in all single cells
    #     variantFile = "dvh-UA3-152-03-singlecell-variants-2callers-80percent-2cells_noan-bed.txt"
    # print(bamheader, fastaFile)

    cmd = "%s %s -l %s -f %s %s > %s" %(program,parameters,variantFile,fastaFile,bamFile, countFile)
    print(cmd)
    print
    os.system(cmd)


def parseCounts():
    # open file for reading
    g = open(countFile, 'r')
    h = open(resultFile, 'w')
    alllines = []
    for row in g:
        alllines.append(row)
    myset = set(alllines)
    # readlines
    for line in myset:
        vector = line.split("\t")
        chrname = vector[0]
        chromosome = chrname
        coordinate = vector[1]
        reference = vector[2]
        mutation = "%s-%s-%s" %(chromosome, coordinate, reference)
        depth = int(vector[3])
        depth1 = depth + 1
        referenceDepth = "%s:%s" %(reference, depth)

        # alternative proportion based
        percentA = round((int(vector[5].split(":")[1])+1) / depth1, 3)
        percentC = round((int(vector[6].split(":")[1])+1) / depth1, 3)
        percentG = round((int(vector[7].split(":")[1])+1) / depth1, 3)
        percentT = round((int(vector[8].split(":")[1])+1) / depth1, 3)

        try:
            percentIN = round((int(vector[10].split(":")[1])+1) / depth1, 3)
        except IndexError:
            percentIN = 0

        print(depth1)
        print("A:%s - C:%s - G:%s - T:%s - INS:%s") %(percentA,percentC,percentG,percentT,percentIN)

        mylist = [percentA, percentC, percentG, percentT, percentIN]
        highest = max(mylist)
        #print(highest)
        if depth == 0:
            topbase = "NA"
            nstatus = "NA"
            proportion = 0
            print("Zer0 depth in this region")
        else:
            if highest == percentA:
                topbase = "A"
                print("Highest is: A")
            elif highest == percentC:
                topbase = "C"
                print("Highest is: C")
            elif highest == percentG:
                topbase = "G"
                print("Highest is: G")
            elif highest == percentT:
                topbase = "T"
                print("Highest is: T")
            else:
                highest == percentIN
                topbase = "INDEL" #vector[10].split(":")[0]
                print("Highest is: INDEL")

            if topbase == reference:
                print("Top hit is reference")
                print("PercentINT is %s") %percentIN
                if percentIN < 0.5:
                    nstatus = 0
                    proportion = "%s" %(highest)
                else:
                    nstatus = 1
                    proportion = "%s" %(highest)
            else:
                nstatus = 1
                proportion = "%s" %(highest)

        # read numbers based
        baseA = vector[5].split(":")[0] + ":" + vector[5].split(":")[1]
        baseC = vector[6].split(":")[0] + ":" + vector[6].split(":")[1]
        baseG = vector[7].split(":")[0] + ":" + vector[7].split(":")[1]
        baseT = vector[8].split(":")[0] + ":" + vector[8].split(":")[1]
        try:
            baseIN = vector[10].split(":")[0] + ":" + vector[10].split(":")[1]
        except IndexError:
            baseIN = "NA:0"


        #base = referenceDepth.split(":")[0]
        if depth <= 3 :
            status = "NA"
            alternativeCount = "NA"
        else:
            for b in [baseA, baseC, baseG, baseT, baseIN]:
                base = b.split(":")[0]
                bdepth = b.split(":")[1]

                if b == referenceDepth:
                    # check if there are insertions/deletions
                    try:
                        rdepth = int(vector[10].split(":")[1])
                    except IndexError:
                        rdepth = "NA"

                    if int(depth)-4 <= rdepth <= int(depth):
                        status = 1
                        ins = vector[10].split(":")[0]
                        alternativeCount = "%s:%s" %(ins,rdepth)
                        break
                    else:
                        status = 0
                        alternativeCount = "%s:%s" %(base,depth)
                    break
                else:
                    if base == reference:
                        try:
                            insdepth = int(vector[10].split(":")[1])
                        except IndexError:
                            insdepth = 0
                        if base == "A":
                            rdepth = int(vector[5].split(":")[1])
                            print(rdepth)
                            print(int(rdepth)-4)
                            if int(depth)-4 <= rdepth <= int(depth):
                                status = 0
                                alternativeCount = "A:%s" %(rdepth)
                                break
                            else:
                                status = 1
                                alternativeCount = "A:%s" %(rdepth)
                        elif base == "C":
                            rdepth = int(vector[6].split(":")[1])
                            if int(depth)-4 <= rdepth <= int(depth):
                                status = 0
                                alternativeCount = "C:%s" %(rdepth)
                                break
                            else:
                                status = 1
                                alternativeCount = "C:%s" %(rdepth)
                        elif base == "G":
                            rdepth = int(vector[7].split(":")[1])
                            if int(depth)-4 <= rdepth <= int(depth):
                                status = 0
                                alternativeCount = "G:%s" %(rdepth)
                                break
                            else:
                                status = 1
                                alternativeCount = "G:%s" %(rdepth)
                        elif base == "T":
                            rdepth = int(vector[8].split(":")[1])
                            if int(depth)-4 <= rdepth <= int(depth):
                                status = 0
                                alternativeCount = "T:%s" %(rdepth)
                                break
                            else:
                                status = 1
                                alternativeCount = "T:%s" %(rdepth)
                        else:
                            base == "NA"
                            rdepth = int(vector[10].split(":")[1])
                            if int(depth)-4 <= rdepth <= int(depth):
                                status = 0
                                alternativeCount = "INS:%s" %(rdepth)
                                break
                            else:
                                status = 1
                                alternativeCount = "INS:%s" %(rdepth)
                    else:
                        obase = [x for x in ["A", "C", "G", "T", "NA"] if x != base]
                        for mbase in obase:
                            if mbase == "A":
                                rdepth = int(vector[5].split(":")[1])
                                if int(depth)-4 <= rdepth <= int(depth):
                                    status = 1
                                    alternativeCount = "A:%s" %(rdepth)
                                    break
                                else:
                                    status = 0
                                    alternativeCount = "A:%s" %(rdepth)
                            elif mbase == "C":
                                rdepth = int(vector[6].split(":")[1])
                                if int(depth)-4 <= rdepth <= int(depth):
                                    status = 1
                                    alternativeCount = "C:%s" %(rdepth)
                                    break
                                else:
                                    status = 0
                                    alternativeCount = "C:%s" %(rdepth)
                            elif mbase == "G":
                                rdepth = int(vector[7].split(":")[1])
                                if int(depth)-4 <= rdepth <= int(depth):
                                    status = 1
                                    alternativeCount = "G:%s" %(rdepth)
                                    break
                                else:
                                    status = 0
                                    alternativeCount = "G:%s" %(rdepth)
                            elif mbase == "T":
                                rdepth = int(vector[8].split(":")[1])
                                if int(depth)-4 <= rdepth <= int(depth):
                                    status = 1
                                    alternativeCount = "T:%s" %(rdepth)
                                    break
                                else:
                                    status = 0
                                    alternativeCount = "T:%s" %(rdepth)
                            else:
                                mbase == "NA"
                                try:
                                    rdepth = int(vector[10].split(":")[1])
                                except IndexError:
                                    rdepth = 0
                                if int(depth)-4 <= rdepth <= int(depth):
                                    status = 1
                                    alternativeCount = "INS:%s" %(rdepth)
                                    break
                                else:
                                    status = 0
                                    alternativeCount = "INS:%s" %(rdepth)
        #print("%s %s %s %s %s %s %s Status: %s") %(chromosome, coordinate, referenceDepth, baseA, baseC, baseG, baseT, status)
        line2write = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tStatus:\t%s\t%s\t%s\ttopbase:%s\tproportion:%s\n" %(chromosome, coordinate, referenceDepth, baseA, baseC, baseG, baseT, baseIN, status, alternativeCount, nstatus, topbase, proportion)
        h.write(line2write)

# loop through each single cell folder to read bmreadcount files and parse
n = 1
totalfolders = len(folders)
for folder in folders:
        print("--------------------%s/%s--------------------") %(n,totalfolders)
        print("Processing %s...") %(folder)
        print
        # get count file
        bamFile = glob.glob(folder + "*_marked.bam")[0]
        print(bamFile)
        sample = bamFile.split("_marked.bam")[0]
        countFile = sample + "_2cells_80percent_noan_bamreadcount.txt"
        resultFile = countFile.split("_2cells_80percent_noan_bamreadcount.txt")[0] + "__2cells_80percent_noan_bamreadcount_parsed.txt"
        print(countFile)
        print(resultFile)
        print
        print("Running bamreadcount on %s") %(folder)
        print
        runBamReadCount()
        print
        print("Parsing bamreadcounts in %s") %(folder)
        print
        parseCounts()
        n = n + 1
        #sys.exit()
