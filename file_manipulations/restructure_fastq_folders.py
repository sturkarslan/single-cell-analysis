import os, glob, sys

datafolder = "/Users/sturkars/Documents/Basespace/20180507_Aotwell-76197121/FASTQ_Generation"
newdatafolder = "/Users/sturkars/Documents/Basespace/20180507_Aotwell-76197121"
folders = glob.glob('%s/*' %datafolder)
#print(folders)
folderNames = []
for folder in folders:
    folderNames.append(folder.split('/')[7].split('_')[0])
folderNames = set(folderNames)
#print(folderNames)

for name in folderNames:
    cmd = 'mkdir %s/%s' %(newdatafolder, name)
    print(cmd)
    os.system(cmd)
    files = glob.glob('%s/*/%s*' %(datafolder,name))
    #print(files)
    for file in files:
        cmd2 = 'mv %s %s/%s' %(file, newdatafolder,name)
        print(cmd2)
        os.system(cmd2)
        
