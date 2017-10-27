# Collect list of bam files and write them into a file
import glob, sys, os, string
resultsDir = 'results-09'
organism = "dvh"
outfile = '%s/%s/bamfileList.txt' %(resultsDir, organism)

bamFiles = glob.glob('%s/*/*/*_sorted.bam' %(resultsDir))
print bamFiles

f = open(outfile, 'w')
for file in bamFiles:
    f.write(file+"\n")
