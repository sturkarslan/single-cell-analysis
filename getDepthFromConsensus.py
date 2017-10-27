from __future__ import division
import glob, sys, os, string

folders = glob.glob("results-03/dvh/DvH*/")
# file to get list of all variants in all single cells
variantFile = "dvh-UA3-152-03-singlecell-variants-2callers-80percent-2cells_norepeats-bed.txt"
# # create a list to append Chromosome-coordinate-allele
# variantList = []
# f = open(variantFile, 'r')
# for line in f:
#     fields = line.split("\t")
#     one = fields[0]
#     two = fields[1]
#     three = fields[3]
#     variantID = "%s-%s-%s" %(one, two, three)
#     variantList.append(variantID)

def runBamReadCount():
    program = "~/github/bam-readcount/bin/bam-readcount"
    parameters = "-q 40 -w 0 -b 13"
    fastaFile = "reference/Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.fasta"
    cmd = "%s %s -l %s -f %s %s > %s" %(program,parameters,variantFile,fastaFile,bamFile, countFile)
    print(cmd)
    print
    os.system(cmd)


def parseCounts():
    # open file for reading
    g = open(countFile, 'r')
    h = open(resultFile, 'w')
    # readlines
    for line in g:
        vector = line.split("\t")
        chromosome = vector[0]
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
        print(depth1)
        print("A:%s - C:%s - G:%s - T:%s") %(percentA,percentC,percentG,percentT)

        mylist = [percentA, percentC, percentG, percentT]
        highest = max(mylist)
        print(highest)

        if highest == percentA:
            topbase = "A"
        elif highest == percentC:
            topbase = "C"
        elif highest == percentG:
            topbase = "G"
        else:
            highest == percentT
            topbase = "T"

        if topbase == reference:
            nstatus = 0
            proportion = "%s" %(highest)
        else:
            nstatus = 1
            proportion = "%s" %(highest)

        # read numbers based
        baseA = vector[5].split(":")[0] + ":" + vector[5].split(":")[1]
        baseC = vector[6].split(":")[0] + ":" + vector[6].split(":")[1]
        baseG = vector[7].split(":")[0] + ":" + vector[7].split(":")[1]
        baseT = vector[8].split(":")[0] + ":" + vector[8].split(":")[1]
        #base = referenceDepth.split(":")[0]
        if depth <= 3 :
            status = "NA"
            alternativeCount = "NA"
        else:
            for b in [baseA, baseC, baseG, baseT]:
                base = b.split(":")[0]
                bdepth = b.split(":")[1]

                if b == referenceDepth:
                    status = 0
                    alternativeCount = "%s:%s" %(base,depth)
                    break
                else:
                    if base == reference:
                        if base == "A":
                            rdepth = int(vector[5].split(":")[1])
                            print(rdepth)
                            print(int(rdepth)-4)
                            if int(depth)-4 <= rdepth <= int(depth):
                                status = 0
                                alternativeCount = "A:%s" %(rdepth)
                                break
                            else:
                                status = 0
                                alternativeCount = "A:%s" %(rdepth)
                        elif base == "C":
                            rdepth = int(vector[6].split(":")[1])
                            if int(depth)-4 <= rdepth <= int(depth):
                                status = 0
                                alternativeCount = "C:%s" %(rdepth)
                                break
                            else:
                                status = 0
                                alternativeCount = "C:%s" %(rdepth)
                        elif base == "G":
                            rdepth = int(vector[7].split(":")[1])
                            if int(depth)-4 <= rdepth <= int(depth):
                                status = 0
                                alternativeCount = "G:%s" %(rdepth)
                                break
                            else:
                                status = 0
                                alternativeCount = "G:%s" %(rdepth)
                        else:
                            base == "T"
                            rdepth = int(vector[8].split(":")[1])
                            if int(depth)-4 <= rdepth <= int(depth):
                                status = 0
                                alternativeCount = "T:%s" %(rdepth)
                                break
                            else:
                                status = 0
                                alternativeCount = "T:%s" %(rdepth)
                    else:
                        obase = [x for x in ["A", "C", "G", "T"] if x != base]
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
                            else:
                                mbase == "T"
                                rdepth = int(vector[8].split(":")[1])
                                if int(depth)-4 <= rdepth <= int(depth):
                                    status = 1
                                    alternativeCount = "T:%s" %(rdepth)
                                    break
                                else:
                                    status = 0
                                    alternativeCount = "T:%s" %(rdepth)
        #print("%s %s %s %s %s %s %s Status: %s") %(chromosome, coordinate, referenceDepth, baseA, baseC, baseG, baseT, status)
        line2write = "%s\t%s\t%s\t%s\t%s\t%s\t%s\tStatus:\t%s\t%s\t%s\ttopbase:%s\tproportion:%s\n" %(chromosome, coordinate, referenceDepth, baseA, baseC, baseG, baseT, status, alternativeCount, nstatus, topbase, proportion)
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
        sample = bamFile.split("_marked.bam")[0]
        countFile = sample + "_2cells_norepeats_bamreadcount.txt"
        resultFile = countFile.split("_2cells_norepeats_bamreadcount.txt")[0] + "_2cells_norepeats_bamreadcount_parsed.txt"
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
