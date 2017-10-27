import glob, sys, os, string
organism = 'dvh'
resultsDir = 'results-03/%s'%organism

outfile = '%s/bamFileList-for-monovar.txt' %(resultsDir)
bamFiles = glob.glob('%s/*/*_marked.bam' %(resultsDir))

g = open(outfile, 'w')

for file in bamFiles:
    g.write(file)
    g.write('\n')
    
    

#
# samtools mpileup -BQ0 -d10000 -f ref.fa -q 40 -b filenames.txt | monovar.py -p 0.002 -a 0.2 -t 0.05 -m 2 -f ref.fa -b filenames.txt -o output.vcf
#
# samtools mpileup -BQ0 -d10000 -f reference/Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.fasta -q 40 -b results-09/dvh/bamFileList-for-monovar.txt | monovar.py -p 0.002 -a 0.2 -t 0.05 -m 2 -f reference/Desulfovibrio_vulgaris_str_hildenborough.GCA_000195755.1.30.dna.genome.fasta -b results-09/dvh/bamFileList-for-monovar.txt -o results-09/dvh/monovar-output.vcf