# /lustre/mib-cri/carrol09/python/PythonInstall/bin/python2.7 /lustre/mib-cri/carrol09/Work/MyPipe/Process10/PythonScripts/ExtractSequences.py /lustre/mib-cri/carrol09/Work/PipelinePracticeSet/20111109_RossAdams_DN_HNF1bChIP/Peaks/Macs_Peaks/SLX-4497.739.s_2.bwa.homo_sapiens_Processed_summits.bed /lustre/mib-cri/carrol09/Work/MyPipe/Genomes/GRCh37/homo_sapiens.fa /lustre/mib-cri/carrol09/Work/PipelinePracticeSet/20111109_RossAdams_DN_HNF1bChIP/OutPutSequences.txt


import pysam
import argparse
import textwrap
import os
import ConfigParser
import sys
import subprocess
import re

Locations = sys.argv[1]
Fasta = sys.argv[2]
OutFile = sys.argv[3]
Size = sys.argv[4] 
#Locations = "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/20111109_RossAdams_DN_HNF1bChIP/Peaks/Macs_Peaks/SLX-4497.739.s_2.bwa.homo_sapiens_Processed_summits.bed"
#Fasta = "/lustre/mib-cri/carrol09/Work/MyPipe/Genomes/GRCh37/homo_sapiens.fa"
#OutFile =  "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/20111109_RossAdams_DN_HNF1bChIP/OutPutSequences.txt"

#print(Locations)

FastaFile = pysam.Fastafile(Fasta)
MyFasta = open(OutFile,'w')
bed = open(Locations,"r")
bedList = []
newSetTemp = []
bedListName = []

for range in bed:
	ChompRange = range.rstrip("\n")
	coords = re.split("\t",ChompRange)
	newSetTemp = []
	newSetTemp.append(coords[0])
	newSetTemp.append(((float(coords[1])+float(coords[2]))/2)-int(Size))
	newSetTemp.append(((float(coords[1])+float(coords[2]))/2)+int(Size))
	newSetTemp.append(coords[3])	
	bedList.append(newSetTemp)

bed.close()

K= 0
Missed = 0
for region in bedList:
	K=K+1
	#print(str(region[0]),int(region[1]),int(region[2]))
	Sequence = FastaFile.fetch(str(region[0]),int(region[1]),int(region[2]))
	if len(Sequence) == (int(Size)*2):
		#print(Sequence)
		MyFasta.write(">"+region[3]+"\n")
		MyFasta.write(Sequence+"\n")
	else:
		 Missed = Missed+1

print("Peaks with no sequence found "+str(Missed)+"\n")
MyFasta.close()
