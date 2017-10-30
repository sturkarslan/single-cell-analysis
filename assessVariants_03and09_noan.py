
############# Serdar Turkarslan / Institute for Systems Biology ###############
# Last update: 10/18/2017
###############################################################################
# This script uses a variant file to check each mutation in that
# variant file against all the variations found in single cells
# mutations are assigned status 1: present, 0: there are enough reads but
# no support for mutations 3: not enough reads for mutation present/absent
# callers Script goes through variant cals in each single cell sequencing
# folders and then check bam readcount files to see if there is support.
###############################################################################


from __future__ import division
import glob, sys

folders = glob.glob("results-0*/dvh/DvH*/")
variantFile = "dvh-UA3-152-03and09-singlecell-variants-2callers-80percent-2cells_noan-bed.txt"
fastaFile = "reference/Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.fasta"
outfile = 'dvh-UA3-152-03and09_noan_mutation_counts_verified.txt'
mutationnames = 'dvh-UA3-152-03and09_2ca80p2cenoan_mutation_names.txt'
cellnames = 'dvh-UA3-152-03and09_2ca80p2cenoan_cell_names.txt'
mutationmatrix = 'dvh-UA3-152-03and09_2ca80p2cenoan_mutation_matrix.txt'

# file to get list of all variants in all single cells
# count files
print("Procesing counts files...\n")
countfiles = []
for folder in folders:
    myfile = glob.glob(folder + "*__2cells_80percent_noan_bamreadcount_parsed.txt")
    countfiles.append(myfile)

# Open variant file and loop through each variant to create a variantList
print("Procesing variants file...\n")
f = open(variantFile, 'r')
variantList = []
for line in f:
    fields = line.split("\t")
    fields[-1] = fields[-1].replace('\n','')
    one = fields[0]
    two = fields[1]
    three = fields[3]
    locus = fields[6]
    names = "%s-%s-%s" %(one, two,locus)
    variantID = "%s-%s-%s_%s" %(one, two, three,locus)
    variantList.append(variantID)

v = open(mutationmatrix, 'w') ## open matrix output only file for writing
h = open(outfile, 'w') ##open output file for matrix and cel and mutation names

# get cell names
headerlist = []
headerlist.append("Variant")
for cfile in countfiles:
    headername = cfile[0].split("/")[3].split("__2cells_80percent_noan_bamreadcount_parsed.txt")[0]
    headerlist.append(headername)
headerlist.append("\n")
header = "\t".join(headerlist)
h.write(header)

## write variants into mutation names file
t = open(mutationnames, 'w')
for mutation in variantList:
    mutationShort = mutation.replace("Chromosome", "")
    mutationShort = mutationShort.replace("pDV", "p")
    mutationShort = mutationShort.replace("DVU_", "DVU")
    mutationsecond = mutationShort.split("_")[1] + "_" + mutationShort.split("-")[1] + mutationShort.split("-")[0]
    t.write(mutationsecond)
    t.write("\n")
t.close()

## write cell names into cell names fileName
s = open(cellnames, 'w')
headerlist.pop(0)
swrite = "\t".join(headerlist) + "\n"
s.write(swrite)
s.close()

# loop through each variant and search each count file to collect mutations status
n = 1
for variant2 in variantList:
    totalvariants = len(variantList)
    print("--------------------%s/%s--------------------" %(n,totalvariants)) 
    print("Processing %s..." %(variant2)) 
    print
    statusList = []
    matrixList = []
    variant = variant2.split("_")[0]
    statusList.append(variant2)

    # loop through each single cell folder to read bmreadcount files and parse
    m = 1
    totalfolders = len(folders)
    for folder in folders:
            print("----------%s/%s----------" %(m,totalfolders)) 
            print("Processing %s..." %(folder)) 
            print
            # get count file
            bamFile = glob.glob(folder + "*_marked.bam")[0]
            consensusVariant = glob.glob(folder + "*consensus.variants.FINAL.txt")[0]
            sample = bamFile.split("_marked.bam")[0]
            countFile = sample + "_2cells_80percent_noan_bamreadcount.txt"
            resultFile = countFile.split("_2cells_80percent_noan_bamreadcount.txt")[0] + "__2cells_80percent_noan_bamreadcount_parsed.txt"

            # create a list of variants from consesnus variant calls for each folder
            consensusList = []
            k = open(consensusVariant, 'r')
            for line in k:
                row = line.split("\t")
                chromosome = row[0]
                coord = row[1]
                reference = row[2]
                joint = "%s-%s-%s" %(chromosome, coord, reference)
                consensusList.append(joint)

            # create a dictionary of mutations as keys and status as values
            countDict = {}
            g = open(resultFile, 'r')
            for line in g:
                row = line.split("\t")
                chromosome = row[0]
                coord = row[1]
                reference = row[2]
                status = row[11]
                joint2 = "%s-%s-%s" %(chromosome, coord, reference)
                joint2 = joint2.split(":")[0]
                countDict[joint2] = status
            keylist = list(countDict.keys())

            # check variants if they are in consensus, if not check in bam counts if not assign 3 (NA)
            if variant in consensusList:
                status = "1"
            elif variant in keylist:
                status = countDict[variant]
            else:
                status = "3"

            m = m + 1
            statusList.append(status)
            matrixList.append(status)
            
    n = n + 1
    line2write = "\t".join(statusList) + "\n" ## for matrix and names
    vwrite = "\t".join(matrixList) + "\n" ## for matrix ony
    h.write(line2write)
    v.write(vwrite)
