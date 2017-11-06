#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 11:27:14 2017

@author: sturkars
"""
from __future__ import division
from subprocess import Popen, PIPE, STDOUT
import glob, sys, os, string, subprocess

runningDir = "/Users/sturkars/github/single-cell-analysis/"

paths = ["/Volumes/omics4tb/sturkarslan/dvh-coculture-rnaseq/dvh-single-cells/results-03/dvh/DvH_03*/",
         "/Volumes/omics4tb/sturkarslan/dvh-coculture-rnaseq/dvh-single-cells/results-09/dvh/DvH_09*/",
         "/Volumes/omics4tb/sturkarslan/EPD/evolved_lines/after_300g/results/dvh/*/",
         "/Volumes/omics4tb/sturkarslan/clonal-isolates/results/dvh/*/"]

folders = []
for path in paths:
    folder = glob.glob(path)
    folders.append(folder)
# create a flat list out of nsted lists
folderlist = [item for sublist in folders for item in sublist]

n = 1
totalfolders = len(folderlist)
for folder in folderlist:
    print("--------------------%s/%s-------------------- %(n,totalfolders)") 
    print("Processing %s... %(folder)") 
    print
    #get count file
    bamFile = glob.glob(folder + "*_marked.bam")[0]
    print(bamFile)
    sample = bamFile.split("_marked.bam")[0]
    countFile = sample + "_allsamples_bamreadcount.txt"
    jobscriptfile = sample + "_jobscript.csh"
    logfile = sample + "_jobscriptlog.txt"
    resultFile = countFile.split("_allsamples_bamreadcount.txt")[0] + "_allsamples_bamreadcount_parsed.txt"
     
    # write to job file
    with open(jobscriptfile,'w') as g:
       g.write('#!/bin/bash\n\n')
       g.write('#$ -N allsamples-%s\n'  %(n))
       g.write('#$ -o %s\n' %(logfile))
       g.write('#$ -e %s\n' %(logfile))
       g.write('#$ -P Bal_sturkars\n')
       g.write('#$ -pe serial 4\n')
       g.write('#$ -q baliga\n')
       g.write('#$ -S /bin/bash\n\n')
    
       g.write('#Sample: %s\n' %(sample))
        
       # changing terminal to bash
       g.write('bash\n\n')
        
       # changing to the appropriate directory
       g.write('cd %s\n\n' %(runningDir)) 
    
       # command for running conversion from tab-delimited to fastq file
       g.write('echo "## Running getDepthFromConsensus_allsamples.py" \n') 
        
       g.write('python getDepthFromConsensus_allsamples.py %s %s %s %s' %(folder,bamFile,countFile,resultFile))            
    
    g.close()
   
    n = n + 1
    sys.exit()
