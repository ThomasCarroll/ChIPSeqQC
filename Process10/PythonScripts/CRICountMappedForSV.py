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
OutFile = sys.argv[2]

TempNames = GetDirGetName(BamFile)
pysam.index(BamFile)
SamFile = pysam.Samfile(BamFile,"rb")
AllMapped = SamFile.mapped
ReAlignLog = open(OutFile,"a")
ReAlignLog.write(str(AllMapped))	
ReAlignLog.write("\n")
ReAlignLog.close()
