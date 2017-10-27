from __future__ import division
import glob, sys, os, string, csv, itertools
import pandas as pd

folders = glob.glob("results-09/dvh/DvH*/")
variantFile = "dvh-UA3-152-09-singlecell-variants-2callers-80percent-2cells_noan-bed.txt"
fastaFile = "reference/Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.fasta"

# file to get list of all variants in all single cells
# count files
countfiles = []
for folder in folders:
    myfile = glob.glob(folder + "/*__2cells_80percent_noan_bamreadcount_parsed.txt")
    countfiles.append(myfile)

#print(countfiles)
#sys.exit()
# outfile
outfile = 'dvh-UA3-152-09_noan_mutation_counts.txt'

# Open variant file and loop through each variant to create a variantList
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
variantList = set(variantList)
#print(variantList)
#sys.exit()

h = open(outfile, 'w')
headerlist = []
headerlist.append("Variant")
#headerlist.append("Name")
for cfile in countfiles:
    headername = cfile[0].split("/")[3].split("__2cells_80percent_noan_bamreadcount_parsed.txt")[0]
    #print(headername)
    headerlist.append(headername)
headerlist.append("\n")
header = "\t".join(headerlist)
h.write(header)

# loop through each variant and search each count file to collect mutations status
for variant2 in variantList:
    statusList = []
    variant = variant2.split("_")[0]
    statusList.append(variant2)

    #print("Variant: %s") %(variant)

    for countfile in countfiles:
        #print(countfile)
        g = open(countfile[0], 'r')

        for line in g:
            row = line.split("\t")
            chromosome = row[0]
            coord = row[1]
            reference = row[2]
            joint = "%s-%s-%s" %(chromosome, coord, reference)
            joint = joint.split(":")[0]
            status = row[11]
            #print(joint)

            if variant == joint:
                #print("Found variant!!!    %s | %s") %(variant, joint)
                status = status
                break
            else:
                status = "3"
        statusList.append(status)
    statusList.append("\n")
    line2write = "\t".join(statusList)
    h.write(line2write)
    #sys.exit()
