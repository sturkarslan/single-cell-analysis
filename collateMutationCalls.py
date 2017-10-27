from __future__ import division
import glob, sys, os, string, csv, itertools
import pandas as pd

folders = glob.glob("results-03/dvh/DvH*/")
# file to get list of all variants in all single cells
variantFile = "dvh-UA3-152-03-singlecell-variants-2callers-80percent-2cells_norepeats-bed.txt"
# count files
countfiles = glob.glob("results-03/dvh/DvH*/*_2cells_norepeats_bamreadcount_parsed.txt")
# outfile
outfile = 'results-03/dvh/dvh-UA3-152-03-mutation_from_bamcounts_2cells_norepeats.txt'

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

h = open(outfile, 'w')
headerlist = []
headerlist.append("Variant")
#headerlist.append("Name")
for cfile in countfiles:
    headername = cfile.split("/")[3].split("_2cells_norepeats_bamreadcount_parsed.txt")[0]
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

    print("Variant: %s") %(variant)

    for countfile in countfiles:
        g = open(countfile, 'r')

        for line in g:
            row = line.split("\t")
            chromosome = row[0]
            coord = row[1]
            reference = row[2]
            joint = "%s-%s-%s" %(chromosome, coord, reference)
            joint = joint.split(":")[0]
            status = row[8]

            if variant == joint:
                #print("Found variant!!!")
                status = status
                break
            else:
                status = "NA"
        statusList.append(status)
    statusList.append("\n")
    line2write = "\t".join(statusList)
    h.write(line2write)
sys.exit()
