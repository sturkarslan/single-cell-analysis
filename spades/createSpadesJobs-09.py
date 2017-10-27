import glob, sys, os, string

################# Main #################    
# Input files    
organism = "dvh"
dataDir = "DvH_single_cell_amplified_genome_strain_UA3_152_09-30575559" 
genomeDir = "reference"
genomeGff = '%s/%s.GCA_000195755.1.30.gtf' %(genomeDir, organism)
resultsDir = "results-09"
fastqcDir = '%s/fastqc' %(resultsDir)
pipelineLog = '%s/pipelineLog.txt' %(resultsDir)
knownSites = '%s-variants-compiled.vcf' %(organism)
topDir = "/proj/omics4tb/sturkarslan/dvh-coculture-rnaseq/dvh-single-cells"

# Process fastQ Files
fastqFoldersAll = glob.glob('%s/*/' %(dataDir))

# remove trimmed directory
fastqFolders=[element for element in fastqFoldersAll if element not in ('DvH_single_cell_amplified_genome_strain_UA3_152_09-30575559/trimmed/')]
print
print 'FASTQ Folders: %s' %len(fastqFolders)

#fastqFolders.remove('DvH_single_cell_amplified_genome_strain_UA3_152_03-30586562/trimmed/')

folderCount = 1
for folder in fastqFolders:
    sampleFolder = folder.split("/")[1] # Folder containing reads
    print '\033[33mProcessing Folder: %s of %s (%s)\033[0m' %(folderCount, len(fastqFolders), sampleFolder )
    fastqFilesFirst = glob.glob('%s/*R1*.fastq.gz' %(folder)) # 1st file the folder
    
    
    for fastqfile in fastqFilesFirst:
        i = fastqfile.split("/")[2]
        print '%s' %(i)
    
    
    filePairs = []
    fileCount = 1
    for file in fastqFilesFirst:
        i = file.split("/")[2]
        print
        print '\033[32m Processing File: %s of %s (%s)\033[0m' %(fileCount, len(fastqFilesFirst), i )
        
        # Collect Sample attributes
        fileName = file.split("_R1")[0]
        sampleResultsDir = resultsDir+ '/'+organism+'/'+fileName.split("/")[1]
        sampleTitle = fileName.split("/")[2] 
        lane = fileName.split("/")[2].split("_")[2]
        sampleName = fileName.split("/")[2]
        #sampleId = fileName.split("/")[2].split("_")[0].split("-")[3] # uncomment for actual run
        sampleId = fileName.split("/")[2].split("_")[0] # uncomment for test run
        # create read pairs
        firstPair = "--pe1-1 " + fileName + "_R1_001.fastq.gz"
        secondPair = "--pe1-2 " + fileName + "_R2_001.fastq.gz"
        combinedPair = firstPair + " " + secondPair
        filePairs.append(combinedPair)
        print "Combined Pair: %s" %(combinedPair)
        
        fileCount = fileCount + 1
    
    # I/O files
    inputFiles = " ".join(filePairs)
    jobFileName = '%s/%s/%s/%s/runSpades.csh' %(topDir, resultsDir, organism, sampleFolder)
    spadesOutputFile = '%s/%s/%s/%s/spades_output' %(topDir, resultsDir, organism, sampleFolder)
    
    # create spades command
    spadesCmd = "/users/sturkars/SPAdes-3.5.0-Linux/bin/spades.py --sc %s -o %s" %(inputFiles, spadesOutputFile)
    print spadesCmd
    
    # write to job file
    with open(jobFileName,'w') as g:
        g.write('#!/bin/bash\n\n')
        g.write('#$ -N spades-%s\n'%(folderCount))
        g.write('#$ -o %s/%s/%s/%s/spades_log.txt\n' %(topDir, resultsDir, organism, sampleFolder))
        g.write('#$ -e %s/%s/%s/%s/spades_log.txt\n' %(topDir, resultsDir, organism, sampleFolder))
        g.write('#$ -P Bal_sturkars\n')
        g.write('#$ -pe serial 4\n')
        g.write('#$ -q baliga\n')
        g.write('#$ -S /bin/bash\n\n')

        g.write('#Sample: %s\n' %(sampleFolder))
        
        # changing terminal to bash
        g.write('bash\n\n')
        
        # Create output directory
        g.write('mkdir %s\n\n' %(spadesOutputFile))

        # changing to the appropriate directory
        g.write('cd %s\n\n' %(topDir)) 

        # command for running conversion from tab-delimited to fastq file
        g.write('echo "## Converting tab-delimited reads to fastq" \n') 
        #g.write('#/tools/bin/python scripts/tagToReadsCollectPlate3.py %s %s\n\n' %(str(vector[1]), outputFileName))
        
        # spades command
        g.write('%s' %(spadesCmd))
        

    g.close()
    folderCount = folderCount + 1
    
    # submit each job with qsub
    cmd = 'qsub %s' %jobFileName
    print cmd
    print
    os.system(cmd)
