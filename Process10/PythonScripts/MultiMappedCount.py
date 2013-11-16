import pysam
import re
import sys
import Pyrex
from datetime import *
import os
from operator import itemgetter

def GetDirGetName(BamFile):
	fileSimple = re.split(".bam",re.split("/.*/",BamFile)[1])[0]
	dir = re.findall("/.*/",BamFile)[0]
	return [dir,fileSimple]


BamFile = sys.argv[1]
for agr in sys.argv:
	print(agr)

TempNames = GetDirGetName(BamFile)
OutFileCentral = str(TempNames[0])+str(TempNames[1])+"_MultiMappedDist.txt"
OutFileCentral2 = str(TempNames[0])+str(TempNames[1])+"_MultiMappedCount.txt"
Count = open(OutFileCentral,"wb")	
Count2 = open(OutFileCentral2,"wb")	

NMapping = {}
print 'Creating Dup dictionary...'
MultiOrNot = 0
#BamFile = "/lustre/mib-cri/carrol09/Work/JoesChristopher/bioinf/MappingChecking/50/FQDir/SampleNewSet_90_50.bwa.RealignedMM9.bam"
#BamFile = sys.argv[1]
print(BamFile)
pysam.index(BamFile)
SamFile = pysam.Samfile(BamFile,"rb")
for alignedread in SamFile.fetch():
        for tag, value in alignedread.tags:
            if tag == 'X0':
                if float(value) > 0:
                	MultiOrNot +=1
                if value not in NMapping.keys():
                	NMapping[value] = 0
                NMapping[value]+=1

          ##
Count.write(str(MultiOrNot))
for key in sorted(map(float,NMapping.keys())):
	Count2.write(str(key)+"\t"+str(NMapping[key]))  	
          

Count.close()
Count2.close()



