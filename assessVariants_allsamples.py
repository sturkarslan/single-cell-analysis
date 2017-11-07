
############# Serdar Turkarslan / Institute for Systems Biology ###############
# Last update: 11/07/2017
###############################################################################
# This script uses a variant file to check each mutation in that
# variant file against all the variations found in single cells
# mutations are assigned status 1: present, 0: there are enough reads but
# no support for mutations 3: not enough reads for mutation present/absent
# callers Script goes through variant cals in each single cell sequencing
# folders and then check bam readcount files to see if there is support.
# this script works on merged files to verify mutations in all samples.
###############################################################################
from __future__ import division
import glob, sys, os, string, csv, itertools
#import pandas as pd

variantFile = "/Volumes/omics4tb/sturkarslan/dvh-coculture-rnaseq/dvh-single-cells/dvh-UA3-152-03and09-singlecell-variants-2callers-80percent-2cells_noan-bed.txt"
fastaFile = "/Volumes/omics4tb/sturkarslan/dvh-coculture-rnaseq/dvh-single-cells/reference/Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.fasta"
outfile = '/Volumes/omics4tb/sturkarslan/dvh-mutation-verifications/dvh-UA3-152-03and09_allsamples_mutation_counts_verified.txt'
mutationnames = '/Volumes/omics4tb/sturkarslan/dvh-mutation-verifications/dvh-UA3-152-03and09_allsamples_mutation_names.txt'
cellnames = '/Volumes/omics4tb/sturkarslan/dvh-mutation-verifications/dvh-UA3-152-03and09_allsamples_cell_names.txt'
mutationmatrix = '/Volumes/omics4tb/sturkarslan/dvh-mutation-verifications/dvh-UA3-152-03and09_allsamples_mutation_matrix.txt'

## paths for all sample folders
paths = ["/Volumes/omics4tb/sturkarslan/dvh-coculture-rnaseq/dvh-single-cells/results-03/dvh/DvH_03*/",
         "/Volumes/omics4tb/sturkarslan/dvh-coculture-rnaseq/dvh-single-cells/results-09/dvh/DvH_09*/",
         "/Volumes/omics4tb/sturkarslan/EPD/evolved_lines/after_300g/results/dvh/*/",
         "/Volumes/omics4tb/sturkarslan/clonal-isolates/results/dvh/*/"]
# create a list of all folders from these paths
folders = []
for path in paths:
    folder = glob.glob(path)
    folders.append(folder)
# create a flat list out of nested lists
folderlist = [item for sublist in folders for item in sublist]


# file to get list of all variants in all single cells
# count files
print("Procesing counts files...\n")
countfiles = []
for folder in folderlist:
    myfile = glob.glob(folder + "*_allsamples_bamreadcount.txt")
    countfiles.append(myfile)
#print(countfiles)

# Open variant file and loop through each variant to create a variantList
print("Procesing variants file...\n")
f = open(variantFile, 'r')
variantList = []
for line in f:
    fields = line.split("\t")
    one = fields[0]
    two = fields[1]
    three = fields[3]
    locus = fields[6]
    names = "%s-%s-%s" %(one, two,locus)
    variantID = "%s-%s-%s_%s" %(one, two, three,locus)
    variantList.append(variantID)


v = open(mutationmatrix, 'w') ## open matrix output only file for writing
h = open(outfile, 'w') ##open output file for matrix and cel and mutation names
headerlist = []
headerlist.append("Variant")
#headerlist.append("Name")
for cfile in countfiles:
    #print(cfile)
    headername = cfile[0].split("/")[3].split("_allsamples_bamreadcount_parsed.txt")[0]
    #print(headername)
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
    mutationShort = mutationShort.split("_")[1] + "." + mutationShort.split("-")[1] + mutationShort.split("-")[0]
    t.write(mutationShort)
    t.write("\n")
#twrite = "\n".join(variantList)
t.close()

## write cell names into cell names fileName
s = open(cellnames, 'w')
swrite = "\t".join(headerlist)
s.write(swrite)
s.close()

n = 1
# loop through each variant and search each count file to collect mutations status
for variant2 in variantList:
    totalvariants = len(variantList)
    print("--------------------%s/%s--------------------") %(n,totalvariants)
    print("Processing %s...") %(variant2)
    print
    statusList = []
    matrixList = []
    variant = variant2.split("_")[0]
    statusList.append(variant2)
    #print(statusList)

    # loop through each single cell folder to read bmreadcount files and parse
    m = 1
    totalfolders = len(folderlist)
    for folder in folderlist:
            print("----------%s/%s----------") %(m,totalfolders)
            print("Processing %s...") %(folder)
            print
            # get count file
            bamFile = glob.glob(folder + "*_marked.bam")[0]
            consensusVariant = glob.glob(folder + "*consensus.variants.FINAL.txt")[0]
            #print(bamFile, variantFile)
            sample = bamFile.split("_marked.bam")[0]
            countFile = sample + "_allsamples_bamreadcount.txt"
            resultFile = countFile.split("_allsamples_bamreadcount.txt")[0] + "_allsamples_bamreadcount_parsed.txt"
            #print(consensusVariant)

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
            sys.exit()
    n = n + 1
    #print(statusList)
    line2write = "\t".join(statusList) + "\n" ## for matrix and names
    vwrite = "\t".join(matrixList) + "\n" ## for matrix ony
    h.write(line2write)
    v.write(vwrite)
