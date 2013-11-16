getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
Args <- commandArgs(trailingOnly = TRUE)
library(raster)
source("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/RScripts/Workflow_Functions.r")


ConfigFile <- readIniFile(file.path(getwd(),"Temp","config.ini"))


LocationsDir <- ConfigFile[ConfigFile[,2] %in% "tempdirectory",3]
WkgDir<- ConfigFile[ConfigFile[,2] %in% "workingdirectory",3]
BamDir<- ConfigFile[ConfigFile[,2] %in% "bamdirectory",3]
FQDir<- ConfigFile[ConfigFile[,2] %in% "fastqdirectory",3]
Genome <- ConfigFile[ConfigFile[,2] %in% "genome",3]
GenomeFileOptions <-  ConfigFile[ConfigFile[,1] %in% "Genomes",]
GenomeFile <- GenomeFileOptions[GenomeFileOptions[,2] %in% tolower(Genome),3]


sampleSheet <- read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=",")
SamplesAndInputs <- matrix(data=c(gsub("_Processed\\.bam","",sampleSheet[,"Processed_bamFileName"]),gsub("_Processed\\.bam","",sampleSheet[match(sampleSheet[,"InputToUse"],sampleSheet[,"SampleName"],incomparable=NA),"Processed_bamFileName"])),ncol=2,byrow=F,dimnames=list(NULL,c("Samples","Inputs")))
SamplesAndInputs <- SamplesAndInputs[!is.na(SamplesAndInputs[,"Inputs"]) & ! as.vector(SamplesAndInputs[,"Samples"]) %in% "No_Processed_Bam",]

sampleSheet <- read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=",")
SamplesAndInputs2 <- matrix(data=c(gsub("_Processed\\.bam","",sampleSheet[,"Processed_bamFileName"]),gsub("_Processed\\.bam","",sampleSheet[match(sampleSheet[,"InputToUse"],sampleSheet[,"GenomicsID"],incomparable=NA),"Processed_bamFileName"])),ncol=2,byrow=F,dimnames=list(NULL,c("Samples","Inputs")))
SamplesAndInputs2 <- SamplesAndInputs2[!is.na(SamplesAndInputs2[,"Inputs"]) & ! as.vector(SamplesAndInputs2[,"Samples"]) %in% "No_Processed_Bam",]

Temp <- rbind(SamplesAndInputs,SamplesAndInputs2)
Temp <- Temp[match(unique(Temp[,1]),Temp[,1]),]

SamplesAndInputs <- matrix(Temp,ncol=2,byrow=F)

Config <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/Config/Config.txt",sep="\t",header=F)

CallPeaks <- ConfigFile[ConfigFile[,2] %in% "callmacspeaks",3]

if(CallPeaks %in% "Yes"){
	CallPeaks <- TRUE
}else{
	CallPeaks <- FALSE	
}

if(nrow(SamplesAndInputs) != 0 & CallPeaks){
MacsGenomes <- matrix(c("HG18","hs","GRCh37","hs","MM9","mm"),ncol=2,byrow=T)
SicerGenomes <- matrix(c("HG18","hg18","GRCh37","hg19","MM9","mm9"),ncol=2,byrow=T)
TPICsGenomes <- matrix(c("HG18","hg18","GRCh37","GRCh37","MM9","mm9"),ncol=2,byrow=T)


MFOLD_PARAMETER <- as.vector(Config[Config[,1] == "Macs_mFold",2])
FRAGMENTLENGTH_PARAMETER <- as.numeric(as.vector(Config[Config[,1] == "Fragment_Length",2]))
SHIFTSIZE_PARAMETER <- round(FRAGMENTLENGTH_PARAMETER/2)
GENOME_PARAMETER <- Genome
MACSGENOME_PARAMETER <- as.vector(MacsGenomes[MacsGenomes[,1] == Genome,2]) 
SICERGENOME_PARAMETER <- as.vector(SicerGenomes[SicerGenomes[,1] == Genome,2]) 
TPICSGENOME_PARAMETER <- as.vector(TPICsGenomes[TPICsGenomes[,1] == Genome,2]) 


if(GENOME_PARAMETER %in% "GRCh37"){
	GENOME_LENGTHS <- "/lustre/mib-cri/carrol09/MyPipe/bedFiles/hg19.txt"
}
if(GENOME_PARAMETER %in% "HG18"){
	GENOME_LENGTHS <- "/lustre/mib-cri/carrol09/MyPipe/bedFiles/hg18.txt"
}


TempFull <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/MultiMacs.xml",header=F,stringsAsFactors=F,quote="")
TempFull[grep("RunDirectory_ToChange",TempFull[,1]),1] <- gsub("RunDirectory_ToChange",getwd(),TempFull[grep("RunDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("BamDirectory_ToChange",TempFull[,1]),1] <- gsub("BamDirectory_ToChange",file.path(getwd(),"bamFiles"),TempFull[grep("BamDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1] <- gsub("WorkingDirectory_ToChange",getwd(),TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("GenomeLengths_ToChange",TempFull[,1]),1] <- gsub("GenomeLengths_ToChange",GENOME_LENGTHS,TempFull[grep("GenomeLengths_ToChange",TempFull[,1]),1])


SpecialisationSet <- TempFull[grep("specialisation identifier",TempFull[,1]):max(grep("specialisation>",TempFull[,1])),1]
AllSets <- rep(list(SpecialisationSet),nrow(SamplesAndInputs))
BigSet <- vector("character")
for(i in 1:length(AllSets)){
  AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])] <- gsub("Indentifier_ToChange",getRandString(len=6),AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("Test_ToChange",AllSets[[i]])] <- gsub("Test_ToChange",paste(SamplesAndInputs[i,1],"_Processed",sep=""),AllSets[[i]][grep("Test_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("Control_ToChange",AllSets[[i]])] <- gsub("Control_ToChange",paste(SamplesAndInputs[i,2],"_Processed",sep=""),AllSets[[i]][grep("Control_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("GM_ToChange",AllSets[[i]])] <- gsub("GM_ToChange",MACSGENOME_PARAMETER,AllSets[[i]][grep("GM_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("GS_ToChange",AllSets[[i]])] <- gsub("GS_ToChange",SICERGENOME_PARAMETER,AllSets[[i]][grep("GS_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("GT_ToChange",AllSets[[i]])] <- gsub("GT_ToChange",TPICSGENOME_PARAMETER,AllSets[[i]][grep("GT_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("MF_ToChange",AllSets[[i]])] <- gsub("MF_ToChange",MFOLD_PARAMETER,AllSets[[i]][grep("MF_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("SS_ToChange",AllSets[[i]])] <- gsub("SS_ToChange",SHIFTSIZE_PARAMETER,AllSets[[i]][grep("SS_ToChange",AllSets[[i]])])
  BigSet <- c(BigSet,AllSets[[i]])
}

FirstAllExceptSpecialisations <- TempFull[c(1:grep("specialisations",TempFull[,1])[1]),1]
SecondAllExceptSpecialisations <- TempFull[c(grep("specialisations",TempFull[,1])[2]:nrow(TempFull)),1]
Tommy <- c(FirstAllExceptSpecialisations,BigSet,SecondAllExceptSpecialisations)
write.table(Tommy,file=file.path(LocationsDir,"TrialMACS.xml"),qmethod="escape",quote=F,sep="\t",col.names=F,row.names=F)
cat("Submitting jobs!............")
system(paste("java -jar /lustre/mib-cri/carrol09/MyPipe/workflow-all-1.2-SNAPSHOT.jar --mode=lsf ",file.path(LocationsDir,"TrialMACS.xml"),sep=""),wait=TRUE,intern=FALSE)
cat("Jobs Submitted!\n")


MacsFiles <- dir(path=file.path(WkgDir,"Peaks","Macs_Peaks"),pattern="*peaks.xls$")
PositivePeaks <- MacsFiles[-grep("negative",MacsFiles)]

#sampleSheet <- read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=",",stringsAsFactors=F)
sampleSheet <- ReadAndLock(file.path(WkgDir,"SampleSheet.csv"),WkdDir,SAF=F,napTime=5)

for(i in 1:nrow(sampleSheet)){
	if(gsub("_Processed.bam","",sampleSheet[i,"Processed_bamFileName"]) %in% gsub("_Processed_peaks.xls","",PositivePeaks)){
		MacsToRead <- file.path(WkgDir,"Peaks","Macs_Peaks",PositivePeaks[gsub("_Processed_peaks.xls","",PositivePeaks) %in% gsub("_Processed.bam","",sampleSheet[i,"Processed_bamFileName"])])
		print(MacsToRead)
		DataIn <- read.delim(MacsToRead,sep="\t",comment.char="#")
		sampleSheet[i,"MacsPeaks"] <- nrow(DataIn)		
		sampleSheet[i,"Macs_name"] <- MacsToRead
		print(nrow(DataIn))
	}
}

MacsGenes <- dir(path=file.path(WkgDir,"Peaks","Macs_Peaks"),pattern="*Annotated.xls$")
for(i in 1:nrow(sampleSheet)){
	if(gsub("_Processed.bam","",sampleSheet[i,"Processed_bamFileName"]) %in% gsub("_Processed_peaks_Annotated.xls","",MacsGenes)){
		MacsToAnnotate <- file.path(WkgDir,"Peaks","Macs_Peaks",MacsGenes [gsub("_Processed_peaks_Annotated.xls","",MacsGenes) %in% gsub("_Processed.bam","",sampleSheet[i,"Processed_bamFileName"])])
		sampleSheet[i,"Macs_Genes_Annotation"] <- MacsToAnnotate
	}
}

MacsGenes <- dir(path=file.path(WkgDir,"Peaks","Macs_Peaks"),pattern="*Annotated.summary$")
for(i in 1:nrow(sampleSheet)){
	if(gsub("_Processed.bam","",sampleSheet[i,"Processed_bamFileName"]) %in% gsub("_Processed_peaks_Annotated.summary","",MacsGenes)){
		MacsToAnnotateSum <- file.path(WkgDir,"Peaks","Macs_Peaks",MacsGenes [gsub("_Processed_peaks_Annotated.summary","",MacsGenes) %in% gsub("_Processed.bam","",sampleSheet[i,"Processed_bamFileName"])])
		DataInAnno <- read.delim(MacsToAnnotateSum,sep="\t",comment.char="#",header=T)[,1]
		sampleSheet[i,"Macs_Genes_In_Peaks"] <- DataInAnno[1]
		sampleSheet[i,"Macs_Peaks_In_Genes"] <- DataInAnno[2]
	}
}
MacsGO <- dir(path=file.path(WkgDir,"Peaks","Macs_Peaks"),pattern="*_GO_Results.txt$")
for(i in 1:nrow(sampleSheet)){
	if(gsub("_Processed.bam","",sampleSheet[i,"Processed_bamFileName"]) %in% gsub("_Processed_GO_Results.txt","",MacsGO)){
		MacsGOToAnnotateSum <- file.path(WkgDir,"Peaks","Macs_Peaks",MacsGO [gsub("_Processed_GO_Results.txt","",MacsGO) %in% gsub("_Processed.bam","",sampleSheet[i,"Processed_bamFileName"])])
		DataInAnno <- read.delim(MacsGOToAnnotateSum,sep="\t",comment.char="#",header=T)
		sampleSheet[i,"Macs_GO_Terms"] <- nrow(DataInAnno[DataInAnno[,"Adjusted_P_Value"] < 0.01 & !is.na(DataInAnno[,"Adjusted_P_Value"]),])
		sampleSheet[i,"Macs_GO_Table"] <- MacsGOToAnnotateSum
	}

}



#write.table(sampleSheet,file=file.path(WkgDir,"SampleSheet.csv"),row.names=F,sep=",")
WriteAndUnlock(sampleSheet,file.path(WkgDir,"SampleSheet.csv"))	
}else{cat("No Peaks To Call")}

write.table("Complete",file.path(LocationsDir,paste(Args[1],"_MainMacsProcess.txt",sep="")),col.names=T,row.names=F,sep=",",quote=F)


