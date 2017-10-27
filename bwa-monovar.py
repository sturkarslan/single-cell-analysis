import glob, sys, os, string

# Process fastQ Files
def runPipeline():
    fastqFolders = glob.glob('%s/*/' %(dataDir))
    print
    print 'FASTQ Folders: %s' %fastqFolders
    print len(fastqFolders)
    # remove trimmed directory
    fastqFolders.remove('DvH_single_cell_amplified_genome_strain_UA3_152_09-30575559/trimmed/')
    print 'New FASTQ Folders: %s' %fastqFolders
    print len(fastqFolders)
    for folder in fastqFolders:
        sampleFolder = folder.split("/")[1] # Folder containing reads
        #print sampleFolder
        fastqFilesFirst = glob.glob('%s/*R1*.fastq.gz' %(folder)) # 1st file the folder
        print
        print 'FASTQ Files: %s' %fastqFilesFirst
        
        for file in fastqFilesFirst:
            print
            print "\033[33m Processing %s \033[0m" %(file)
            print
            fileName = file.split("_R1")[0]
            print fileName
            
            # Collect Sample attributes
            sampleResultsDir = resultsDir+ '/'+organism+'/'+fileName.split("/")[1]
            sampleTitle = fileName.split("/")[1] 
            lane = fileName.split("/")[2].split("_")[1]
            sampleName = fileName.split("/")[2]
            #sampleId = fileName.split("/")[2].split("_")[0].split("-")[3] # uncomment for actual run
            sampleId = fileName.split("/")[2].split("_")[0] # uncomment for test run
            print
            print "FileName: %s Lane: %s sample: %s ID: %s" %(fileName, lane, sampleName, sampleId)
            
            firstPair = fileName + "_R1_001.fastq.gz"
            secondPair = fileName + "_R2_001.fastq.gz"
            print
            print "First Pair: %s, Second Pair: %s" %(firstPair, secondPair)
            
            # Run pipeline commands
            # Run Fastqc
            runQC(firstPair, secondPair)
            
            # Run trimmomatic
            firstPairTrimmedPaired, secondPairTrimmedPaired = runTrim(firstPair, secondPair)
            
            # Run bwa    
            runBWA(firstPairTrimmedPaired,  secondPairTrimmedPaired, lane, sampleName, sampleId, sampleResultsDir, sampleTitle)
            
            # Run samtools fixmate
            runSamtoolsFixmate(sampleResultsDir,sampleTitle)
            
            # Run samtools sort
            runSamtoolsSort(sampleResultsDir,sampleTitle)
            
            #sys.exit()
            
        
    return firstPair, secondPair, libraryName, sampleResultsDir


# Quality control
def runQC(firstPair, secondPair):
    print
    print "\033[34m Running FastQC Quality Control \033[0m"
    
    # create results folder
    if not os.path.exists('%s' %(fastqcDir)):
        print '\033[31m %s directory does NOT exists.Creating one! \033[0m' %(fastqcDir)  
        os.makedirs('%s' %(fastqcDir))
    else:
        print '\033[31m %s directory exists. Not creating. \033[0m' %(fastqcDir)
    # run Command and write output into both screen and logfile with 2>&1 | tee -a %s
    cmd = '%s -t 4 -o %s %s %s' %(fastqc, fastqcDir, firstPair, secondPair)
    print 'FastQC Command:', cmd
    #os.system(cmd)
    print
    print
#
#
#
#
# Trim reads
def runTrim(firstPair, secondPair):
    print
    print "\033[34m Running Read Trimming... \033[0m"
    # Program Parameters
    illuminaClip = "ILLUMINACLIP:/users/sturkars/Trimmomatic-0.35/adapters/NexteraPE-PE.fa:2:30:10" #Remove adapters 
    leading = "LEADING:3" #Remove leading low quality or N bases
    trailing = "TRAILING:3" #Remove trailing low quality or N bases
    slidingWindow = "SLIDINGWINDOW:4:20" #Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 20
    minLen = "MINLEN:36"
    paired = "PE" # or SE for single end
    threads = "-phred33 -threads 8"
      
    # define result files
    filesFolder = firstPair.split('/')[0]
    firstPairTrimmedPaired = filesFolder+"/trimmed/"+firstPair.split('/')[2].split(".fastq.gz")[0] + "_paired_trimmed.fastq.gz"
    secondPairTrimmedPaired = filesFolder+"/trimmed/"+secondPair.split('/')[2].split(".fastq.gz")[0] + "_paired_trimmed.fastq.gz"
    firstPairTrimmedUnpaired = filesFolder+"/trimmed/"+firstPair.split('/')[2].split(".fastq.gz")[0] + "_unpaired_trimmed.fastq.gz"
    secondPairTrimmedUnpaired = filesFolder+"/trimmed/"+secondPair.split('/')[2].split(".fastq.gz")[0] + "_unpaired_trimmed.fastq.gz"
    trimDir = filesFolder+"/trimmed/"
    
    # create trim folder
    if not os.path.exists('%s' %(trimDir)):
        print '\033[31m %s directory does NOT exists.Creating one! \033[0m' %(trimDir)  
        os.makedirs('%s' %(trimDir))
    else:
        print '\033[31m %s directory exists. Not creating. \033[0m' %(trimDir)
    
    # define command
    cmd = '/users/sturkars/java/bin/java -Xmx128m -jar %s %s %s %s %s %s %s %s %s %s %s %s %s %s 2>&1 | tee -a %s' %(trimmomaticPath, paired, threads, firstPair, secondPair, firstPairTrimmedPaired, firstPairTrimmedUnpaired, secondPairTrimmedPaired, secondPairTrimmedUnpaired, illuminaClip, leading, trailing, slidingWindow, minLen, pipelineLog)
    print "Trimmomatic Command: ", cmd
    #os.system(cmd)
    return firstPairTrimmedPaired, secondPairTrimmedPaired
    print
    print



    
def runBWA(firstPairTrimmedPaired, secondPairTrimmedPaired,lane,sampleName,sampleId,sampleResultsDir,sampleTitle):
    print
    print "\033[34m Running BWA alignment... \033[0m"
    
    # create results folder
    if not os.path.exists('%s' %(sampleResultsDir)):
        print '\033[31m %s directory does NOT exists.Creating one! \033[0m' %(sampleResultsDir)  
        os.makedirs('%s' %(sampleResultsDir))
    else:
        print '\033[31m %s directory exists. Not creating. \033[0m' %(sampleResultsDir)
    
    # modify read group information
    readGroup = "@RG\\tID:%s\\tPL:ILLUMINA\\tSM:%s\\tLB:%s" %(sampleId, sampleName, lane)
    # bwa run command
    cmd = "%s mem -R '%s' %s %s %s > %s/%s.sam" %(bwaPath, readGroup, genomeFasta, firstPairTrimmedPaired, secondPairTrimmedPaired, sampleResultsDir, sampleTitle)
    print "Run BWA Command", cmd
    #os.system(cmd)
    print
 
  
  
      
def runSamtoolsFixmate(sampleResultsDir,sampleTitle):
    print
    print "\033[34m Running SAMtools fixmate... \033[0m"
    # fixmate run command
    cmd = '%s fixmate -O bam %s/%s.sam %s/%s_fixmate.bam 2>&1 | tee -a %s' %(samtoolsPath, sampleResultsDir, sampleTitle, sampleResultsDir, sampleTitle, pipelineLog)
    print "Samtools Fixmate Command: ", cmd
    #os.system(cmd)



    
def runSamtoolsSort(sampleResultsDir,sampleTitle):
    print
    print "\033[34m Running SAMtools sort.. \033[0m"
    # run command
    cmd = '%s sort -O bam -o %s/%s_sorted.bam -T tmp/%s_temp %s/%s_fixmate.bam 2>&1 | tee -a %s' %(samtoolsPath, sampleResultsDir, sampleTitle, sampleTitle, sampleResultsDir, sampleTitle, pipelineLog)
    print "Samtools Sort Command: ", cmd
    # index bam file
    cmd2 = '%s index %s/%s_sorted.bam 2>&1 | tee -a %s' %(samtoolsPath, sampleResultsDir, sampleTitle,  pipelineLog)
    #os.system(cmd)
    #os.system(cmd2)






def genomeIndexes():
    print
    print "\033[34m Creating Genome indexes... \033[0m"
    # indexing genomes
    cmd1=bwaPath+' index '+genomeFasta
    cmd2=samtoolsPath+' faidx '+genomeFasta
    cmd3='%s -Xmx4G -jar %s CreateSequenceDictionary R=%s O=%s.dict 2>&1 | tee -a %s' %(javaPath,piccardPath,genomeFasta,genomeFasta, pipelineLog)
    print "BWA Genome index command: ", cmd1
    print
    print "Samtools Genome index Command: ", cmd2
    print
    print "GATK CreateDictionary Command: ", cmd3

    #os.system(cmd1)
    #os.system(cmd2)
    #os.system(cmd3)

    
################# Main #################    
# Input files    
organism = "dvh"
dataDir = "DvH_single_cell_amplified_genome_strain_UA3_152_09-30575559" 
genomeDir = "reference"
genomeFasta = '%s/Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.fasta' %(genomeDir)
genomeGff = '%s/%s.GCA_000195755.1.30.gtf' %(genomeDir, organism)
resultsDir = "results-09"
fastqcDir = '%s/fastqc' %(resultsDir)
pipelineLog = '%s/pipelineLog.txt' %(resultsDir)
knownSites = '%s-variants-compiled.vcf' %(organism)
# snpEff databases
if organism == "mmp":
    snpEffDatabase = "Methanococcus_maripaludis_S2_uid58035"
if organism == "dvh":
    snpEffDatabase = "dvh-genome"    

### Create sequnce dictioonary
#/usr/bin/java -Xmx4G -jar /users/sturkars/picard-tools-1.139/picard.jar CreateSequenceDictionary REFERENCE=Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.fasta OUTPUT=Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.fasta.dict

# Programs
bwaPath = "/users/sturkars/bwa-0.7.12/bwa" # path to STAR executable
samtoolsPath = "/users/sturkars/samtools-1.2/bin/samtools"   # path to Trimmomatic executable
bcftoolsPath = "/users/sturkars/bcftools-1.1/bcftools"
fastqc = "/users/sturkars/FastQC/fastqc" # path to fastqc executable
javaPath = "/usr/bin/java"
piccardPath = "/users/sturkars/picard-tools-1.139/picard.jar"
trimmomaticPath = "/users/sturkars/Trimmomatic-0.35/trimmomatic-0.35.jar"
gatkPath = "/users/sturkars/gatk/GenomeAnalysisTK.jar"
tabixPath = "/users/sturkars/htslib-1.1/tabix"
varscanPath = "/users/sturkars/VarScan.v2.3.9.jar"
snpEffPath = "/users/sturkars/snpEff/snpEff.jar"

# Run functions
genomeIndexes()
firstPair, secondPair, libraryName, sampleResultsDir = runPipeline()
    