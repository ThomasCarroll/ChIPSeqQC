import re
import sys
import Pyrex
from datetime import *
import os
from operator import itemgetter


BedFile = sys.argv[1]
BedFile2 = sys.argv[2]

bed = open(BedFile,"r")
bed2 = open(BedFile2,"w")

count=1
for line in bed:
	ChompRange = line.rstrip("\n")
	coords = re.split("\t",ChompRange)
	bed2.write(str(coords[0])+"\t"+str(coords[1])+"\t"+str(coords[2])+"\t"+"Peak_"+str(count))
	#print(len(coords))
	#print(range(3,int(len(coords))))
	BedLineLength = int(len(coords))
	for n in range(3,BedLineLength):
		bed2.write("\t"+str(coords[n]))	
	bed2.write("\n")
	count = count+1
bed.close()
bed2.close()




