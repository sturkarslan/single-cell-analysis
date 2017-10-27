import glob, sys, os, string

folders = glob.glob("results-03/dvh/DvH*/")
outfile = "dvh-UA3-152-03_bamfiles.txt"
f = open(outfile, "w")

bamfileList = []
for folder in folders:
    bamFile = glob.glob(folder + "*_marked.bam")[0]
    bamfileList.append(bamFile)

line2write = "\n".join(bamfileList)
f.write(line2write)    
