# /lustre/mib-cri/carrol09/python/PythonInstall/bin/python2.7 /lustre/mib-cri/carrol09/Work/MyPipe/Process10/PythonScripts/ExtractSequences.py /lustre/mib-cri/carrol09/Work/PipelinePracticeSet/20111109_RossAdams_DN_HNF1bChIP/Peaks/Macs_Peaks/SLX-4497.739.s_2.bwa.homo_sapiens_Processed_summits.bed /lustre/mib-cri/carrol09/Work/MyPipe/Genomes/GRCh37/homo_sapiens.fa /lustre/mib-cri/carrol09/Work/PipelinePracticeSet/20111109_RossAdams_DN_HNF1bChIP/OutPutSequences.fa


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

FastaFile = pysam.Fastafile(Fasta)
MyFasta = open(OutFile,'w')
bed = open(Locations,"r")
bedList = []
newSetTemp = []
bedListName = []



K= 0
Missed = 0
for region in bed:
	ChompRange = region.rstrip("\n")
	coords = re.split("\t",ChompRange)
	print(str(coords[0]),int(coords[1]),int(coords[2]))
	Sequence = FastaFile.fetch(str(coords[0]),int(coords[1]),int(coords[2]))
	MyFasta.write(">"+coords[4]+"\n")
	MyFasta.write(Sequence+"\n")
MyFasta.close()
bed.close()
