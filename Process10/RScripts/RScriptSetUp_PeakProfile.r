getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))

Args <- commandArgs(trailingOnly = TRUE)
library(raster)
library(XML)
#ConfigFile <- readIniFile(file.path(Args[2],"config.ini"))
ConfigFile <- readIniFile(file.path(getwd(),"Temp","config.ini"))
LocationsDir <- ConfigFile[ConfigFile[,2] %in% "tempdirectory",3]
WkgDir<- ConfigFile[ConfigFile[,2] %in% "workingdirectory",3]
BamDir<- ConfigFile[ConfigFile[,2] %in% "bamdirectory",3]
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


CallPeaks <- ConfigFile[ConfigFile[,2] %in% "callmacspeaks",3]
if(CallPeaks %in% "Yes"){
	CallPeaks <- TRUE
}else{
	CallPeaks <- FALSE	
}

if(nrow(SamplesAndInputs) != 0 & CallPeaks){


GENOME_PARAMETER <- ConfigFile[ConfigFile[,2] %in% "genome",3] 
#Config <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/Config/Config.txt",sep="\t",header=F)

#GENOME_PARAMETER <- as.vector(Config[Config[,1] == "Genome",2])

if(tolower(GENOME_PARAMETER) %in% tolower("GRCh37")){
	GENOME_LENGTHS <- "/lustre/mib-cri/carrol09/MyPipe/bedFiles/hg19.txt"
}
if(tolower(GENOME_PARAMETER) %in% tolower("HG18")){
	GENOME_LENGTHS <- "/lustre/mib-cri/carrol09/MyPipe/bedFiles/hg18.txt"
}

sampleSheet <- read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=",")
MacsPeaks <- sampleSheet[!is.na(sampleSheet[,"Macs_name"]),"Macs_name"]
SampleNames <- sampleSheet[!is.na(sampleSheet[,"Macs_name"]),"GenomicsID"]
bamNames <- sampleSheet[!is.na(sampleSheet[,"Macs_name"]),"Processed_bamFileName"]
MacsPeaks <- gsub(".xls",".bed",MacsPeaks)




TempFull <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/PeakProfileMeta.xml",header=F,stringsAsFactors=F,quote="")
TempFull[grep("RunDirectory_ToChange",TempFull[,1]),1] <- gsub("RunDirectory_ToChange",getwd(),TempFull[grep("RunDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("BamDirectory_ToChange",TempFull[,1]),1] <- gsub("BamDirectory_ToChange",BamDir,TempFull[grep("BamDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1] <- gsub("WorkingDirectory_ToChange",WkgDir,TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("FastaFile_ToChange",TempFull[,1]),1] <- gsub("FastaFile_ToChange",GenomeFile,TempFull[grep("FastaFile_ToChange",TempFull[,1]),1])
TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1] <- gsub("TempDirectory_ToChange",LocationsDir,TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("genomeFile_ToChange",TempFull[,1]),1] <- gsub("genomeFile_ToChange",GENOME_LENGTHS,TempFull[grep("genomeFile_ToChange",TempFull[,1]),1])


SpecialisationSet <- TempFull[grep("specialisation identifier",TempFull[,1]):max(grep("specialisation>",TempFull[,1])),1]
AllSets <- rep(list(SpecialisationSet),length(MacsPeaks))
BigSet <- vector("character")
for(i in 1:length(AllSets)){
  AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])] <- gsub("Indentifier_ToChange",getRandString(len=6),AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("Test_ToChange",AllSets[[i]])] <- gsub("Test_ToChange",SampleNames[i],AllSets[[i]][grep("Test_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("PeakFileLocation_ToChange",AllSets[[i]])] <- gsub("PeakFileLocation_ToChange",MacsPeaks[i],AllSets[[i]][grep("PeakFileLocation_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("BamFile_ToChange",AllSets[[i]])] <- gsub("BamFile_ToChange",bamNames[i],AllSets[[i]][grep("BamFile_ToChange",AllSets[[i]])])
  BigSet <- c(BigSet,AllSets[[i]])
}

FirstAllExceptSpecialisations <- TempFull[c(1:grep("specialisations",TempFull[,1])[1]),1]
SecondAllExceptSpecialisations <- TempFull[c(grep("specialisations",TempFull[,1])[2]:nrow(TempFull)),1]
Tommy <- c(FirstAllExceptSpecialisations,BigSet,SecondAllExceptSpecialisations)
write.table(Tommy,file=file.path(LocationsDir,"RunPeakProfiler.xml"),qmethod="escape",quote=F,sep="\t",col.names=F,row.names=F)
cat("Submitting jobs!............")
system(paste("java -Xmx1G -jar /lustre/mib-cri/carrol09/MyPipe/workflow-all-1.2-SNAPSHOT.jar --mode=lsf ",file.path(LocationsDir,"RunPeakProfiler.xml"),sep=""),wait=TRUE,intern=FALSE)
}
write.table("Complete",file.path(LocationsDir,paste(Args[1],"_MainPeakProfileProcess.txt",sep="")),col.names=T,row.names=F,sep=",",quote=F)


