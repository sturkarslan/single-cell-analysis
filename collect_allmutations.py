
############# Serdar Turkarslan / Institute for Systems Biology ###############
# Last update: 11/07/2017
###############################################################################
# Collect all mutations in EPDs and pouplations into a single file
###############################################################################
from __future__ import division
import glob, sys
#import pandas as pd

outfile = '/Volumes/omics4tb/sturkarslan/dvh-mutation-verifications/allmutations_exceptsinglecell.txt'

# paths for all sample folders
paths = [#"/Volumes/omics4tb/sturkarslan/EPD/evolved_lines/after_300g/results/dvh/*/",
         "/Volumes/omics4tb/sturkarslan/EPD/EPD_seq/results/dvh/*/"]
         #"/Volumes/omics4tb/sturkarslan/clonal-isolates/results/dvh/*/"]

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


h = open(outfile, 'w')
#countfiles = []
for folder in folderlist:
    variantfile = glob.glob(folder + "*.consensus.variants.FINAL.txt")
    samplename = variantfile[0].split("/")[8]
    f = open(variantfile[0], 'r')
    # loop through each line
    for line in f:
        fields = line.split("\t")
        chromosome = fields[0]
        coordinate = fields[1]
        alternative = fields[4]
        programs = fields[13]
        frequencies = fields[14]
        variantname = "%s-%s-%s" %(chromosome,coordinate,alternative)
        programcount = len(programs.split(":"))
        print(programs, programcount)
        line2write = variantname + "\t" + samplename + "\t" + line
        h.write(line2write)
h.close()    
    
#
## Open variant file and loop through each variant to create a variantList
#print("Procesing variants file...\n")
#f = open(variantFile, 'r')
#variantList = []
#for line in f:
#    fields = line.split("\t")
#    one = fields[0]
#    two = fields[1]
#    three = fields[3]
#    locus = fields[6]
#    names = "%s-%s-%s" %(one, two,locus)
#    variantID = "%s-%s-%s_%s" %(one, two, three,locus)
#    variantList.append(variantID)
#
#
#v = open(mutationmatrix, 'w') ## open matrix output only file for writing
#h = open(outfile, 'w') ##open output file for matrix and cel and mutation names
#headerlist = []
#headerlist.append("Variant")
#
#for cfile in countfiles:
#    if cfile[0].split("/")[6] == "dvh":
#        headername = cfile[0].split("/")[7].split("_allsamples_bamreadcount_parsed.txt")[0] + "-" + cfile[0].split("/")[4].split("_allsamples_bamreadcount_parsed.txt")[0]
#        
#    if cfile[0].split("/")[7] == "dvh":
#        headername = cfile[0].split("/")[8].split("_allsamples_bamreadcount_parsed.txt")[0].split("-")[0]
#     
#    if cfile[0].split("/")[6] == "after_300g":
#        headername = cfile[0].split("/")[9].split("_allsamples_bamreadcount_parsed.txt")[0].split("_")[0]
#    
#    #print(headername)
#    headerlist.append(headername)
#headerlist.append("\n")
#header = "\t".join(headerlist)
#h.write(header)
#
### write variants into mutation names file
#t = open(mutationnames, 'w')
#for mutation in variantList:
#    mutationShort = mutation.replace("Chromosome", "")
#    mutationShort = mutationShort.replace("pDV", "p")
#    mutationShort = mutationShort.replace("DVU_", "DVU")
#    mutationShort = mutationShort.split("_")[1] + "." + mutationShort.split("-")[1] + mutationShort.split("-")[0]
#    t.write(mutationShort)
#    t.write("\n")
##twrite = "\n".join(variantList)
#t.close()
#
### write cell names into cell names fileName
#s = open(cellnames, 'w')
#swrite = "\t".join(headerlist)
#s.write(swrite)
#s.close()
#
#n = 1
## loop through each variant and search each count file to collect mutations status
#for variant2 in variantList:
#    totalvariants = len(variantList)
#    print("--------------------%s/%s--------------------" %(n,totalvariants)) 
#    print("Processing %s..." %(variant2)) 
#    print
#    statusList = []
#    matrixList = []
#    variant = variant2.split("_")[0]
#    variant2 = variant2.replace("\n", "")
#    statusList.append(variant2)
#    
#    # loop through each single cell folder to read bmreadcount files and parse
#    m = 1
#    totalfolders = len(folderlist)
#    for folder in folderlist:
#            print("----------%s/%s----------" %(m,totalfolders)) 
#            print("Processing %s..." %(folder)) 
#            print
#            # get count file
#            bamFile = glob.glob(folder + "*_marked.bam")[0]
#            consensusVariant = glob.glob(folder + "*consensus.variants.FINAL.txt")[0]
#            sample = bamFile.split("_marked.bam")[0]
#            countFile = sample + "_allsamples_bamreadcount.txt"
#            resultFile = countFile.split("_allsamples_bamreadcount.txt")[0] + "_allsamples_bamreadcount_parsed.txt"
#            
#            # create a list of variants from consesnus variant calls for each folder
#            consensusList = []
#            k = open(consensusVariant, 'r')
#            for line in k:
#                row = line.split("\t")
#                chromosome1 = row[0]
#                if chromosome1 == "NC_002937":
#                    chromosome = "Chromosome"
#                elif chromosome1 == "NC_005863":
#                    chromosome = "pDV"
#                else:
#                    chromosome = chromosome1
#                coord = row[1]
#                reference = row[2]
#                #alt = row[15]
#                joint = "%s-%s-%s" %(chromosome, coord, reference)
#                consensusList.append(joint)
#            
#
#            # create a dictionary of mutations as keys and status as values
#            countDict = {}
#            g = open(resultFile, 'r')
#            for line in g:
#                row1 = line.split("\t")
#                chromosome1 = row1[0]
#                if chromosome1 == "NC_002937":
#                    chromosome = "Chromosome"
#                elif chromosome1 == "NC_005863":
#                    chromosome = "pDV"
#                else:
#                    chromosome = chromosome1
#                coord = row1[1]
#                reference = row1[2]
#                alt = row1[10]
#                status = row1[11]
#                ins = row1[7]
#                joint2 = "%s-%s-%s" %(chromosome, coord, reference)
#                rcount = str(reference.split(":")[1])
#                joint2 = joint2.split(":")[0]
#                summary = reference + "|" + alt + "|" + ins + "|" + rcount + "|" + status
#                countDict[joint2] = summary
#            keylist = list(countDict.keys())
#            #print(keylist)
#
#            # check variants if they are in consensus, if not check in bam counts if not assign 3 (NA)
#            print(variant)
#            ## mutations is called by callers
#            if variant in consensusList:
#                status = countDict[variant] + "|1"
#                print("%s in consensusList" %variant)
#            
#            # mutation is not    
#            elif variant in keylist:
#                print("%s in keyList" %variant)
#                verify = int(countDict[variant].split("|")[3])
#                if verify >= 5:
#                    print("%s is higher than 5" %verify)
#                    status = countDict[variant] + "|0"
#                else:
#                    status = countDict[variant] + "|3"
#                    print("%s is LOWER than 5" %verify)
#                
#            # mutation can not be called
#            else:
#                status = "3"
#            
#            ## Do mutation presence/absence calls agree?
#            if status == "3":
#                status = "3"
#            else:
#                verify1 = status.split("|")[4]
#                verify2 = status.split("|")[5]
#            if verify1 == verify2:
#                status = status + "|OK"
#            else:
#                status = status + "|CHECK"
#            
#            m = m + 1
#            statusList.append(status)
#            matrixList.append(status)
#    #sys.exit()
#    n = n + 1
#    #print(statusList)
#    #statusList.append("\n")
#    #print(statusList[0:3])
#    line2write = "\t".join(statusList)
#    line2write = line2write + "\n"
#    #print(line2write)
#    h.write(line2write)
#    ## for matrix and names
#    vwrite = "\t".join(matrixList) + "\n" ## for matrix ony
#  
#    v.write(vwrite)
#h.close()
#v.close()
