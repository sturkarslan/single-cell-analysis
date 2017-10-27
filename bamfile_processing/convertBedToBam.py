import glob, sys, os, string
resultsDir = 'results-09'
# Collect list of bam files
bamFiles = glob.glob('%s/*/*/*_marked.bam' %(resultsDir))

for bamFile in bamFiles:
    print "Analyzing %s" %bamFile
    bedFile = bamFile.split('_marked.bam')[0] + ".bed"
    cmd1 = "bamToBed -i %s > %s" %(bamFile, bedFile)
    cmd2 = "cp %s %s/bedFiles/" %(bedFile, resultsDir)

    print "Now running %s" %cmd1
    os.system(cmd1)
    print "Now running %s" %cmd2
    os.system(cmd2)
    sys.exit()
    

