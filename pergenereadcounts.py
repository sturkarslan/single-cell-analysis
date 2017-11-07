#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 11:23:44 2017

@author: sturkars
"""
import glob, sys, os

# -----------------//----------------------
gtfFile = "/proj/omics4tb/sturkarslan/dvh-coculture-rnaseq/dvh-single-cells/reference/Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.gtf"

resultsFile = "/proj/omics4tb/sturkarslan/dvh-mutation-verifications/pergenecoverage_all.txt"
    
## paths for all sample folders
paths = ["/proj/omics4tb/sturkarslan/dvh-coculture-rnaseq/dvh-single-cells/results-03/dvh/DvH_03*/",
         "/proj/omics4tb/sturkarslan/dvh-coculture-rnaseq/dvh-single-cells/results-09/dvh/DvH_09*/"]

# create a list of all folders from these paths
folders = []
for path in paths:
    folder = glob.glob(path)
    folders.append(folder)
# create a flat list out of nested lists
folderlist = [item for sublist in folders for item in sublist]
bamfilelist = glob.glob("/proj/omics4tb/sturkarslan/dvh-coculture-rnaseq/dvh-single-cells/results-0*/dvh/*/*_marked.bam")
bamfiles = " ".join(bamfilelist)
cmd = "bedtools multicov -bams %s -bed %s > %s" %(bamfiles,gtfFile,resultsFile)
print(cmd)
os.system(cmd)

