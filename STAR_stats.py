# Collect stats from all bam files and write them into a file
import glob, sys, os, string
resultsDir = '/proj/omics4tb/sturkarslan/icalvum_assembly/results_c5/
outfile = '/proj/omics4tb/sturkarslan/icalvum_assembly/results_c5/collected_STAR_stats.txt'
print('Outfile is: %s' %outfile)
#sys.exit

logFiles = glob.glob('%s/*/*_star_Log.final.out' %(resultsDir))
print(logFiles)
sys.exit()

g = open(outfile, 'w')
#header = 'Organism\tSample\tTotal\tSecondary\tSupplementary\tDuplicates\tMapped\tPaired\tRead1\tRead2\tProperlyPaired\twithItself\tsingletons\tMateMappedDiffChr\tMateMappedDiffChrMapq5\n'
#g.write(header)
for file in logFiles:
    sampleName = file.split('_star_Log.final.out')[-2].split('/')[-1]
    g.write('%s\t'%organism)
    g.write('%s\t'%sampleName)
    f = open(file, 'r')
    lineNo = 0
    for line in f:
        vector = line.split("\n")
        vectorJoined = '\t'.join(vector)
        if lineNo == 4:
            vectorSplit = vectorJoined.split('(')[1].split(":")[0]
            print(vectorSplit)
        else:
            vectorSplit = vectorJoined.split(' ')[0]

        g.write('%s\t'%vectorSplit)
        lineNo = lineNo+1
    g.write("\n")

    #sys.exit()
