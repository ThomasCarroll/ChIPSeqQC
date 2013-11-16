#/lustre/mib-cri/carrol09/python/PythonInstall/bin/python2.7
#OutDir = "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/"
#flag1 = 1

import pysam
import re
import sys
import Pyrex
from datetime import *
import os
from operator import itemgetter
import pysam

def GetDirGetName(BamFile):
	fileSimple = re.split(".bam",re.split("/.*/",BamFile)[1])[0]
	dir = re.findall("/.*/",BamFile)[0]
	return [dir,fileSimple]


BamFile = sys.argv[1] 
OutDir = sys.argv[2]
mapQFilter = int(sys.argv[3])

#if str(sys.argv[4]) == "TRUE":
#	DupFilt = True
#else:
#	DupFilt = False

#BamFile = "/lustre/mib-cri/cclab_stg139/bwa/STG139_N8/STG139_N8.bam"
#OutDir = "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/"
#mapQFilter = 15
#DupFilt = True

SamFile = pysam.Samfile(BamFile,"rb")
TempNames = GetDirGetName(BamFile)

print BamFile
print OutDir
print mapQFilter
#print DupFilt

ChromosomesOfInterest = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chr23","chr24","chrX","chrY","chrM"]
#for k in ChromosomesOfInterest:
#	header3 = SamFile.header
#	outfileName = str(OutDir)+str(TempNames[1])+str(k)+".bam"
#	outfile = pysam.Samfile(outfileName,"wb",header=header3)
#	flagTemp = 0
#	for alignedread in SamFile.fetch(str(k)):
#		if alignedread.cigar is not None:
#			for cigar in alignedread.cigar:
#				if cigar[0] == 4:
#					flagTemp = 1
#			if flagTemp == 1:
#				if flag1 == 1:
#					if alignedread.mapq > 15:
#						outfile.write(alignedread)
#				else:
#					outfile.write(alignedread)
#			flagTemp = 0
#	outfile.close()
#	pysam.index(outfileName)


for k in ChromosomesOfInterest:
	header3 = SamFile.header
	outfileName = str(OutDir)+str(TempNames[1])+"."+str(k)+".bam"
	outfile = pysam.Samfile(outfileName,"wb",header=header3)
	for alignedread in SamFile.fetch(str(k)):
		if alignedread.mapq > mapQFilter and alignedread.mapq is not None:
			if alignedread.is_duplicate == False:
				count = 1
				outfile.write(alignedread)
	outfile.close()
	pysam.index(outfileName)




		

