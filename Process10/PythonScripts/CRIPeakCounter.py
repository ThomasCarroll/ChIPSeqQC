import pysam
import re
import sys
import Pyrex
from datetime import *
import os
from operator import itemgetter

def GetGenome(BamFile):
	genome = "HG18"
	TempHead = GetHeaderInfo(BamFile)
	TempComments = TempHead["Comments"]
	for comments in TempComments:
		if comments.split(":")[0] == "GenomeVersion":
			genome = comments.split(":")[1]
	return(genome)

def GetChromoList():
	#This will need to non-hardcoded and also assembly/species specific
	ChromosomesOfInterest = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chr23","chr24","chrX","chrY","chrM"]
	return ChromosomesOfInterest

def CountRandom(BamFile):
	samIdxStats = pysam.idxstats(BamFile)
	samfile = pysam.Samfile(BamFile,"rb")
	TotalMapped = samfile.mapped
	samfile.close()
	countAlign = 0
	List = GetChromoList()
	for stat in samIdxStats:
		if stat.split()[0] in List:
			MappedforChromosome = stat.split()[2]
			countAlign = countAlign+long(MappedforChromosome)
	RandAlign = TotalMapped-countAlign
	return [BamFile,{"Random":RandAlign}]


def GetDirGetName(BamFile):
	fileSimple = re.split(".bam",re.split("/.*/",BamFile)[1])[0]
	dir = re.findall("/.*/",BamFile)[0]
	return [dir,fileSimple]

def CountPeakList(bedList,filename,outfileName):
	bedfile = open(outfileName2,"w")
	header3 = TempSamFile.header
	outfile = pysam.Samfile(outfileName,"wb",header=header3)
	take = (0, 2, 3) # CIGAR operation (M/match, D/del, N/ref_skip)
	for region in bedList:
		print(str(region[0])+"\t"+str(region[1])+"\t"+str(region[2])+"\n")
		for alignedread in TempSamFile.fetch(str(region[0]),int(region[1]),int(region[2])):
			if printCount == 1000000:
				print(countQC+NotQC)
				printCount = 0	
			if alignedread.mapq > 15:	
				countQC += 1
				if alignedread.is_duplicate == False:
					NonDupCount +=1	
					outfile.write(alignedread)
        				if alignedread.is_reverse: strand = "-"
        				else: strand = "+"
					t = sum([ l for op,l in alignedread.cigar if op in take ])
					Chromosome = Chromosomes[alignedread.tid]
					Start = alignedread.pos
					end = alignedread.pos+t
					dummy = alignedread.mapq
					name = alignedread.qname
					strand = strand
					bedfile.write(str(Chromosome)+"\t"+str(Start)+"\t"+str(end)+"\t"+str(name)+"\t"+str(dummy)+"\t"+str(strand)+"\n")
				else:
					DupCount +=1
			else:
				NotQC += 1
			printCount += 1
	outfile.close()
	TempSamFile.close()
	bedfile.close()
	return([countQC,NotQC,DupCount,NonDupCount])

def GetReadLength(BamFile):
	TempSamFile = pysam.Samfile(BamFile,"rb")
	lineN = 0
	K = []
	for alignedread in TempSamFile.fetch():
		K = alignedread.qlen
		if lineN == 1:
			break
		lineN +=1
	return(K)


def ReformatRegions(BamFile):
	bed = open(BedFile,"r")
	bedList = []
	bedListName = []
	for range in bed:
		ChompRange = range.rstrip("\n")
		coords = re.split("\t",ChompRange)
		newSetTemp = coords
		newSetTemp.append(float(re.sub("chr","",coords[0])))
		newSetTemp[1] = float(newSetTemp[1])+(K+1)
		newSetTemp[2] = float(newSetTemp[2])-(K+1)
		bedList.append(newSetTemp)
		bed.close()
	return(bedList)

def GetHeaderInfo(BamFile):
	samfile = pysam.Samfile(BamFile,"rb")
	header = samfile.header
	HeaderTags = header.keys()
	if "HD" in HeaderTags:
		HDs = header["HD"]
		HDSub = HDs.keys()
		if "SO" in HDSub:
			SortOrder = HDs["SO"]
		else:
			SortOrder = "unknown"
		if "VN"	in HDSub:
			BamVersionNumber = HDs["VN"]	
		else:
			BamVersionNumber = "1.0"
	else:
		SortOrder = "unknown"
		BamVersionNumber = "1.0"
	if "PG" in HeaderTags:
		PGs = header["PG"]
	else:
		PGs = [{'ID': 'Unknown Aligner', 'VN': 'NA'}]
	if "SQ" in HeaderTags:
		SQs = header["SQ"]
	if "RG" in HeaderTags:
		RGs = header["RG"]
	else:
		RGs = []
	if "CO" in HeaderTags:
		COs = header["CO"]
	else:
		COs = []
	HeaderInfo = {"SortOrder":SortOrder,"SamVersion":BamVersionNumber,"ProgramInfo":PGs,"Chromosomes":SQs,"ReadGroups":RGs,"Comments":COs}
	return HeaderInfo



def FixHeaderNew(BamFile):
	OldHeaderInfo = GetHeaderInfo(BamFile)
	samfile = pysam.Samfile(BamFile)
	UpDateRG = True	
	NewHD = {"SO":OldHeaderInfo["SortOrder"],"VN":OldHeaderInfo["SamVersion"]}
	#NewSQDic = {"SQ":dict(zip(samfile.references,samfile.lengths))}
	k = 0
	LongSet = []
	ChromToCheck = GetChromoList()
	for i in samfile.references:
		if samfile.references[k] in ChromToCheck:
			TempDic = {"SN":samfile.references[k],"LN":int(samfile.lengths[k])}
			LongSet.append(TempDic)
		k = k+1		
	NewSQs = LongSet
	NewPGs = OldHeaderInfo["ProgramInfo"]
	NewPGs.append({"ID":"ChIPseqPipeLine","VN":"0.01"})
	NewRGs = []	
	if len(re.findall("SLX",BamFile)) == 1:
		RGID = re.search("SLX-[0-9]*",BamFile).group()
	else:
		RGID = re.split(".bam",re.split("/.*/",BamFile)[1])[0]
	for RGtemp in OldHeaderInfo["ReadGroups"]:
		if "ID" in RGtemp.keys():
			if RGtemp["ID"] == RGID:
				UpDateRG = False
	if UpDateRG == True:
		NewRGs = OldHeaderInfo["ReadGroups"]
		NewRGs.append({"ID":RGID,"SM":RGID})
	else:
		NewRGs = OldHeaderInfo["ReadGroups"]
	NewCO =  OldHeaderInfo["Comments"]
	NewHeader = {"HD":NewHD,"SQ":NewSQs,"RG":NewRGs,"PG":NewPGs,"CO":NewCO}
	return(NewHeader)



BamFile = sys.argv[1]
for agr in sys.argv:
	print(agr)

print(BamFile)

tstart = datetime.now()
SamFile = pysam.Samfile(BamFile,"rb")
TempNames = GetDirGetName(BamFile)
if AllMapped > 0:
	bedList = ReformatRegions(BedFile)
	OutFileCentral = str(TempNames[0])+str(TempNames[1])+"_PeakCounts.txt"
	Metrics = CountPeakList(bedList,BamFile,OutFileCentral,OutFileCentral2)

tend = datetime.now()
TimeStamp = str(date.today())
TimeTaken = str(tend-tstart)
TimeString = TimeStamp+" ["+TimeTaken+"]"
fileLog = open(TempNames[0]+TempNames[1]+"_fileLog.log","wb")
fileLog.write("Genome"+"\t"+"Read Length"+"\t"+"Mapped"+"\t"+"NonRandomChr"+"\t"+"IncludedRegions"+"\t"+"QC > 15"+"\t"+"Unique"+"\t"+"Date Processes [Time processing]"+"\n"+str(genome)+"\t"+str(readlength)+"\t"+str(AllMapped)+"\t"+str(NonRand)+"\t"+str(TotalAfterExclude)+"\t"+str(Metrics[0])+"\t"+str(Metrics[3])+"\t"+TimeString+"\n")
#fileLog.write(str(AllMapped)+"\t"+str(NonRand)+"\t"+str(TotalAfterExclude)+"\t"+str(Metrics[0])+"\t"+str(Metrics[3])+"\n")
fileLog.close()
