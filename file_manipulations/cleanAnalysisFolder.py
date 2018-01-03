##############################################################################
## This script collects all the extra intermediary files from variant calling 
## pipeline for deletion/cleaning up
## last modified: 12/01/2017 by Serdar Turkarslan
##############################################################################

import glob, sys, os, string

resultsDir = '/Volumes/omics4tb/sturkarslan/syntroph-raw-sequences-100_generation/'
folders = glob.glob('%s/*/*' %(resultsDir))

bamFiles = glob.glob('%s/*/*/*.bam*' %(resultsDir))
vcfFiles = glob.glob('%s/*/*/*.vcf' %(resultsDir))
txtFiles = glob.glob('%s/*/*/*.txt' %(resultsDir))
fixmateFiles = glob.glob('%s/*/*/*_fixmate.bam*' %(resultsDir))
recalFiles = glob.glob('%s/*/*/*_recal*' %(resultsDir))
realignFiles = glob.glob('%s/*/*/*_realigned*' %(resultsDir))
intervalFiles = glob.glob('%s/*/*/*.intervals*' %(resultsDir))

#print 'BAM Files: %s' %bamFiles
#print 'VCF Files: %s' %vcfFiles
#print 'TXT Files: %s' %txtFiles

for bam in bamFiles:
    cmd = 'rm %s' %bam
    #print cmd
    #os.system(cmd)
    print

for vcf in vcfFiles:
    cmd = 'rm %s' %vcf
    #print cmd
    #os.system(cmd)
    print

for txt in txtFiles:
    cmd = 'rm %s' %txt
    #print cmd
    #os.system(cmd)
    print

for folder in folders:
    folderName = folder.split('/')[2]
    samName = '%s/%s.sam' %(folder, folderName)
    cmd = 'rm %s' %samName
    #os.system(cmd)

# for fix in fixmateFiles:
#     cmd = 'rm %s' %fix
#     print cmd
#     #os.system(cmd)
#     print

for recal in recalFiles:
    cmd = 'rm %s' %recal
    #print cmd
    #os.system(cmd)
    print

for realign in realignFiles:
    cmd = 'rm %s' %realign
    #print cmd
    #os.system(cmd)
    print

for interval in intervalFiles:
    cmd = 'rm %s' %interval
    print cmd
    os.system(cmd)
    print
