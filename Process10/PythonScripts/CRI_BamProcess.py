import pysam
import re
import sys
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
		if comments.split(":")[0] == "Genome":
			genome = comments.split(":")[1]
	return(genome)

def GetChromoList():
	#This will need to non-hardcoded and also assembly/species specific
	ChromosomesOfInterest = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chr23","chr24","chrX","chrY","chrM"]
	return ChromosomesOfInterest

def CountRandom(BamFile,SequenceDic):
	samIdxStats = pysam.idxstats(BamFile)
	samfile = pysam.Samfile(BamFile,"rb")
	TotalMapped = samfile.mapped
	samfile.close()
	countAlign = 0
	List = ExtractSeqOrder(SequenceDic)
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

def GetNewReferenceOrder(Header):
	List = []
	for ref in HeaderFile["SQ"]:
		List.append(ref["SN"])
	return List

def CountBedList(bedList,filename,outfileName,outfileName2,header4):
	countQC=0
	NotQC=0
	NonDupCount=0
	DupCount=0
	printCount = 0
	NewOrder = GetNewReferenceOrder(header4)
	TempSamFile = pysam.Samfile(filename,"rb")
	Chromosomes = TempSamFile.references
	bedfile = open(outfileName2,"w")
	header3 = TempSamFile.header
	outfile = pysam.Samfile(outfileName,"wb",header=header4)
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
					alignedread.rname = NewOrder.index(Chromosomes[alignedread.rname])
					outfile.write(alignedread)
        				if alignedread.is_reverse: strand = "-"
        				else: strand = "+"
					t = sum([ l for op,l in alignedread.cigar if op in take ])
					Chromosome = NewOrder[alignedread.rname]
					Start = alignedread.pos
					end = alignedread.pos+t
					dummy = alignedread.mapq
					name = alignedread.qname
					strand = strand
					bedfile.write(str(Chromosome)+"\t"+str(Start)+"\t"+str(end)+"\t"+str(name)+"\t"+str(dummy)+"\t"+str(strand)+"\n")
				else:
					DupCount +=1
					if DupFlag == "False":
						alignedread.rname = NewOrder.index(Chromosomes[alignedread.rname])
						outfile.write(alignedread)
        					if alignedread.is_reverse: strand = "-"
        					else: strand = "+"
						t = sum([ l for op,l in alignedread.cigar if op in take ])
						Chromosome = NewOrder[alignedread.rname]
						Start = alignedread.pos
						end = alignedread.pos+t
						dummy = alignedread.mapq
						name = alignedread.qname
						strand = strand
						bedfile.write(str(Chromosome)+"\t"+str(Start)+"\t"+str(end)+"\t"+str(name)+"\t"+str(dummy)+"\t"+str(strand)+"\n")
			else:
				NotQC += 1
				if MapQFlag == "False":	
					NotQC -= 1
					countQC += 1
					if alignedread.is_duplicate == False:
						NonDupCount +=1
						alignedread.rname = NewOrder.index(Chromosomes[alignedread.rname])
						outfile.write(alignedread)
        					if alignedread.is_reverse: strand = "-"
        					else: strand = "+"
						t = sum([ l for op,l in alignedread.cigar if op in take ])
						Chromosome = NewOrder[alignedread.rname]
						Start = alignedread.pos
						end = alignedread.pos+t
						dummy = alignedread.mapq
						name = alignedread.qname
						strand = strand
						bedfile.write(str(Chromosome)+"\t"+str(Start)+"\t"+str(end)+"\t"+str(name)+"\t"+str(dummy)+"\t"+str(strand)+"\n")
					else:
						DupCount +=1
						if DupFlag == "False":
							alignedread.rname = NewOrder.index(Chromosomes[alignedread.rname])
							outfile.write(alignedread)
        						if alignedread.is_reverse: strand = "-"
        						else: strand = "+"
							t = sum([ l for op,l in alignedread.cigar if op in take ])
							Chromosome = NewOrder[alignedread.rname]
							Start = alignedread.pos
							end = alignedread.pos+t
							dummy = alignedread.mapq
							name = alignedread.qname
							strand = strand
							bedfile.write(str(Chromosome)+"\t"+str(Start)+"\t"+str(end)+"\t"+str(name)+"\t"+str(dummy)+"\t"+str(strand)+"\n")
				#else:
				#	countQC += 1
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


def ReformatExcludedRegions(BamFile,Genome,BedFile,SequenceDic):
	K = GetReadLength(BamFile)
	genome = Genome
	ChrOrder = ExtractSeqOrder(SequenceDic)
	bedList = []
	newSetTemp = []
	bedListName = []
	for Chr in ChrOrder:
		print Chr
		bed = open(BedFile,"r")
		for range in bed:
			ChompRange = range.rstrip("\n")
			coords = re.split("\t",ChompRange)
			print(coords[0])
			if coords[0] == Chr:
				print "match"
				newSetTemp = coords
				newSetTemp.append(coords[0])
				newSetTemp[1] = float(newSetTemp[1])+(K+1)
				newSetTemp[2] = float(newSetTemp[2])-(K+1)
				if float(newSetTemp[2]) > float(newSetTemp[1]):
					bedList.append(newSetTemp)
				else:
					print(newSetTemp)
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

def ExtractSeqOrder(SequenceDic):
	SeqsInHeader = pysam.Samfile(SequenceDic,"r")
	chromosomeOrderList = []
	header = SeqsInHeader.header
	SQs = header["SQ"]
	for SQ in SQs:
		chromosomeOrderList.append(SQ["SN"])	
	return chromosomeOrderList 	


def bedListFromSeqDict(BamFile,SequenceDic):
	SeqsInHeader = pysam.Samfile(SequenceDic,"r")
	chromosomeOrderList = []
	header = SeqsInHeader.header
	SQs = header["SQ"]
	for SQ in SQs:
		chromosomeOrderList.append([SQ["SN"],1,SQ["LN"],SQ["SN"]])
	return chromosomeOrderList 		

def FixHeaderNew(BamFile,SequenceDic):
	OldHeaderInfo = GetHeaderInfo(BamFile)
	samfile = pysam.Samfile(BamFile)
	samfileForHeaderOrder = pysam.Samfile(SequenceDic)
	ReferencesFromFile = samfileForHeaderOrder.references
	LengthsFromFile =  samfileForHeaderOrder.lengths
	UpDateRG = True	
	NewHD = {"SO":OldHeaderInfo["SortOrder"],"VN":OldHeaderInfo["SamVersion"]}
	#NewSQDic = {"SQ":dict(zip(samfile.references,samfile.lengths))}
	k = 0
	LongSet = []
	ChromToCheck = ExtractSeqOrder(SequenceDic)
	for i in ReferencesFromFile:
		if ReferencesFromFile[k] in ChromToCheck:
			TempDic = {"SN":ReferencesFromFile[k],"LN":int(LengthsFromFile[k])}
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
#for agr in sys.argv:
#	print(agr)
#
print(BamFile)
#BamFile = "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/20100826_RossInnesC_JC_TumorChIPs/bamFiles/SLX-3230.521.s_6.bwa.homo_sapiens.bam"

#BamFile = "/lustre/mib-cri/carrol09/Work/MyPipe/20111109_RossAdams_DN_HNF1bChIP/OldToTest.bam"
#BamFile = "/lustre/mib-cri/carrol09/Work/MyPipe/20111109_RossAdams_DN_HNF1bChIP/ForTest.bam"
Genome = sys.argv[2]
BedFile = sys.argv[3]
SequenceDic = sys.argv[4]
ExcludedFlag = sys.argv[5]
DupFlag = sys.argv[6]
MapQFlag = sys.argv[7]


#ExcludedFlag = True
#DupFlag = True
#MapQFlag = True



#Genome = "GRCh37"
#BedFile = "/lustre/mib-cri/carrol09/MyPipe/bedFiles/HG19_ExcludedGenome.bed"
#BamFile = "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/Trial/Hickey/20110624_TVAna_CC_MBD24HR/bamFiles/SLX-5380.878.s_1.bwa.homo_sapiens.bam"
#SequenceDic = "/lustre/mib-cri/carrol09/Work/MyPipe/Genomes/GRCh37/GRCh37_ReducedSeqDict.sam"

if ExcludedFlag == "False":
	BedFile = "No_Excluded"



tstart = datetime.now()
pysam.index(BamFile)
SamFile = pysam.Samfile(BamFile,"rb")
AllMapped = SamFile.mapped
Header = SamFile.header
genome = GetGenome(BamFile)
TempNames = GetDirGetName(BamFile)
if AllMapped > 0:
	if BedFile == "No_Excluded":
		bedList = bedListFromSeqDict(BamFile,SequenceDic)
	else:
		bedList = ReformatExcludedRegions(BamFile,Genome,BedFile,SequenceDic)
	
	SamFile.close()
	RandomCount = CountRandom(BamFile,SequenceDic)
	NonRand = AllMapped-RandomCount[1]["Random"]
	OutFileCentral = str(TempNames[0])+str(TempNames[1])+"_Processed.bam"
	OutFileCentral2 = str(TempNames[0])+str(TempNames[1])+"_Processed.bed"
	HeaderFile = FixHeaderNew(BamFile,SequenceDic)
	Metrics = CountBedList(bedList,BamFile,OutFileCentral,OutFileCentral2,HeaderFile)
	TotalAfterExclude = Metrics[0]+Metrics[1]
	readlength = GetReadLength(BamFile)
	pysam.index(OutFileCentral)

if AllMapped == 0:
	readlength = 0
	AllMapped = 0
	NonRand = 0
	TotalAfterExclude = 0
	Metrics = [0,0,0,0]
	OutFileCentral = str(TempNames[0])+str(TempNames[1])+"_Processed.bam"
	OutFileCentral2 = str(TempNames[0])+str(TempNames[1])+"_Processed.bed"
	outfile = pysam.Samfile(OutFileCentral, "wb", header = Header)
	outfile.close()
	pysam.index(OutFileCentral)
	bedfile = open(OutFileCentral2,"w")
	bedfile.write(" "+"\n")
	bedfile.close()

tend = datetime.now()
TimeStamp = str(date.today())
TimeTaken = str(tend-tstart)
TimeString = TimeStamp+" ["+TimeTaken+"]"
fileLog = open(TempNames[0]+TempNames[1]+"_fileLog.log","wb")
fileLog.write("Genome"+"\t"+"Read Length"+"\t"+"Mapped"+"\t"+"NonRandomChr"+"\t"+"IncludedRegions"+"\t"+"QC > 15"+"\t"+"Unique"+"\t"+"Date Processes [Time processing]"+"\n"+str(genome)+"\t"+str(readlength)+"\t"+str(AllMapped)+"\t"+str(NonRand)+"\t"+str(TotalAfterExclude)+"\t"+str(Metrics[0])+"\t"+str(Metrics[3])+"\t"+TimeString+"\n")
#fileLog.write(str(AllMapped)+"\t"+str(NonRand)+"\t"+str(TotalAfterExclude)+"\t"+str(Metrics[0])+"\t"+str(Metrics[3])+"\n")
fileLog.close()
