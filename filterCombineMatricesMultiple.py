
############# Serdar Turkarslan / Institute for Systems Biology ###############
# Last update: 11/01/2017
###############################################################################
# This script filters mutation matrix for mutations that occur in n-numbuer of
# cells. 
###############################################################################

from __future__ import division
import glob


# load variant files
variantfiles = glob.glob("/proj/omics4tb/sturkarslan/dvh-mutation-verifications/dvh-UA3-152-*_curated.txt")
cellfilters = [1,2,5]
nasfilter = [50]

for cellfilter in cellfilters:

    for variantfile in variantfiles:
        samplename = variantfile.split("_")[0]
        cellsfile = samplename + "_SC_curated_" + str(cellfilter) + "cells.txt"

        mutationnames = samplename + "_SC_curated_" + str(cellfilter) + "cells_mutation_names.txt"

        cellnames = samplename + "_SC_curated_" + str(cellfilter) + "cells_cell_names.txt"

        matrixfile = samplename + "_SC_curated_" + str(cellfilter) + "cells_mutation_matrix.txt"

        # Open variant file and loop through each variant to create a variantList
        print("Procesing variants file %s...\n"  %(variantfile))
        print
        f = open(variantfile, 'r')
        headers = f.readline()

        # open output files for writing and add headers
        v = open(cellsfile, 'w')
        #h = open(tencellsfile, 'w')
        k = open(matrixfile, 'w')
        v.write(headers)
        #h.write(headers)

        # start lists to hold 5and 10 cells filtered mutations
        mutations = []
        presents = []
        absents =  []
        nas = []

        # read through each line of file after header andcount presence/absen/na of mutations
        for line in f:
            fields = line.split("\t")
            fields[-1] = fields[-1].replace('\n','')
            matrixonly = fields[1:len(fields)]
            #matrixonly = matrixonly.pop(0)
            print(fields)
            print(matrixonly)
            name = fields[0]

            # count 1s in the line
            present = sum(e.count("1") for e in matrixonly)
            absent = sum(e.count("0") for e in matrixonly)
            na = sum(e.count("3") for e in matrixonly)
            print("Name: %s, Present: %s, Absent: %s, NA: %s \n" %(name, present, absent, na)) 
            presents.append(present)
            absents.append(absent)
            nas.append(na)

            # filter for mutations that are in at least 5 cells
            if(present >= cellfilter):
                mutations.append(name)
                fieldsString = "\t".join(fields)
                print(fieldsString)
                v.write("{}\n".format(fieldsString))
                matrixString = "\t".join(matrixonly)
                k.write("{}\n".format(matrixString))
             

        # ## write variants into mutation names file
        t = open(mutationnames, 'w')
        for mutation in mutations:
            mutationShort = mutation.replace("Chromosome", "")
            mutationShort = mutationShort.replace("pDV", "p")
            mutationShort = mutationShort.replace("DVU_", "DVU")
            mutationShort = mutationShort.split("_")[1] + "_" + mutationShort.split("-")[1] + mutationShort.split("-")[0]
            #mutationShort = mutationShort.split("_")[1] + "_" + mutationShort.split("-")[1] + mutationShort.split("-")[0]
            t.write(mutationShort)
            t.write("\n")
        t.close()

        # ## write cell names into cell names fileName
        s = open(cellnames, 'w')
        cellheader = headers.split("\t")
        cellheader.pop(0)
        swrite = "\t".join(cellheader) + "\n"
        s.write(swrite)
        s.close()
        k.close()
#
#        plt.hist(nas)
#        plt.ylabel('Cells')
#        plt.title('Present')
#        plt.show()

