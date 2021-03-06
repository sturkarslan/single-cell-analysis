# Collect stats from all bam files and write them into a file
import glob, sys, os, string
organism = 'mmp'
resultsDir = '/Volumes/omics4tb/sturkarslan/clonal-isolates/results/%s'%organism
outfile = '/Volumes/omics4tb/sturkarslan/clonal-isolates/results/%s/collected_bam_stats.txt'%organism
print('Outfile is: %s' %outfile)
#sys.exit

statFiles = glob.glob('%s/*/*.flagstats.txt' %(resultsDir))
print(statFiles)

g = open(outfile, 'w')
header = 'Organism\tSample\tTotal\tSecondary\tSupplementary\tDuplicates\tMapped\tPaired\tRead1\tRead2\tProperlyPaired\twithItself\tsingletons\tMateMappedDiffChr\tMateMappedDiffChrMapq5\n'
g.write(header)
for file in statFiles:
    sampleName = file.split('.flagstats.txt')[-2].split('/')[-1]
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
