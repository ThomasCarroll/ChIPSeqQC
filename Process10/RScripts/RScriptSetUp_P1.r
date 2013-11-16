getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))


Args <- commandArgs(trailingOnly = TRUE)
library(raster)
ConfigFile <- readIniFile(file.path(Args[2],"config.ini"))
ConfigFile <- readIniFile(file.path(getwd(),"Temp","config.ini"))
LocationsDir <- ConfigFile[ConfigFile[,2] %in% "tempdirectory",3]
WkgDir<- ConfigFile[ConfigFile[,2] %in% "workingdirectory",3]
BamDir<- ConfigFile[ConfigFile[,2] %in% "bamdirectory",3]
 




#read.delim("/lustre/mib-cri/carrol09/Work/20111109_RossAdams_DN_HNF1bChIP")
sampleSheet <- read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=",")
Samples <- sampleSheet[!is.na(sampleSheet[,"bamFileName"]) & ! sampleSheet[,"Location"] %in% "No_Lims_Location_Known","bamFileName"]
Samples <- gsub("_Processed","",Samples)
Config <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/Config/Config.txt",sep="\t",header=F)

MacsGenomes <- matrix(c("HG18","hs","GRCh37","hs","MM9","mm9"),ncol=2,byrow=T)
SicerGenomes <- matrix(c("HG18","hg18","GRCh37","hg19","MM9","mm9"),ncol=2,byrow=T)
TPICsGenomes <- matrix(c("HG18","hg18","GRCh37","GRCh37","MM9","mm9"),ncol=2,byrow=T)


MFOLD_PARAMETER <- as.vector(Config[Config[,1] == "Macs_mFold",2])
FRAGMENTLENGTH_PARAMETER <- as.numeric(as.vector(Config[Config[,1] == "Fragment_Length",2]))
SHIFTSIZE_PARAMETER <- round(FRAGMENTLENGTH_PARAMETER/2)
GENOME_PARAMETER <- as.vector(Config[Config[,1] == "Genome",2])
MACSGENOME_PARAMETER <- as.vector(MacsGenomes[MacsGenomes[,1] == GENOME_PARAMETER,2]) 
SICERGENOME_PARAMETER <- as.vector(SicerGenomes[SicerGenomes[,1] == GENOME_PARAMETER,2]) 
TPICSGENOME_PARAMETER <- as.vector(TPICsGenomes[TPICsGenomes[,1] == GENOME_PARAMETER,2]) 


if(GENOME_PARAMETER %in% "GRCh37"){
	GENOME_LENGTHS <- "/lustre/mib-cri/carrol09/MyPipe/bedFiles/hg19.txt"
}
if(GENOME_PARAMETER %in% "HG18"){
	GENOME_LENGTHS <- "/lustre/mib-cri/carrol09/MyPipe/bedFiles/hg18.txt"
}


TempFull <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/MultiBamMeta_P1.xml",header=F,stringsAsFactors=F,quote="")
TempFull[grep("RunDirectory_ToChange",TempFull[,1]),1] <- gsub("RunDirectory_ToChange",getwd(),TempFull[grep("RunDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("BamDirectory_ToChange",TempFull[,1]),1] <- gsub("BamDirectory_ToChange",BamDir,TempFull[grep("BamDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1] <- gsub("WorkingDirectory_ToChange",WkgDir,TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("GenomeLengths_ToChange",TempFull[,1]),1] <- gsub("GenomeLengths_ToChange",GENOME_LENGTHS,TempFull[grep("GenomeLengths_ToChange",TempFull[,1]),1])
TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1] <- gsub("TempDirectory_ToChange",LocationsDir,TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1])

SpecialisationSet <- TempFull[grep("specialisation identifier",TempFull[,1]):max(grep("specialisation>",TempFull[,1])),1]
AllSets <- rep(list(SpecialisationSet),length(Samples))
BigSet <- vector("character")
for(i in 1:length(AllSets)){
  AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])] <- gsub("Indentifier_ToChange",getRandString(len=6),AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("Test_ToChange",AllSets[[i]])] <- gsub("Test_ToChange",paste(gsub(".bam","",Samples[i]),sep=""),AllSets[[i]][grep("Test_ToChange",AllSets[[i]])])
  BigSet <- c(BigSet,AllSets[[i]])
}

FirstAllExceptSpecialisations <- TempFull[c(1:grep("specialisations",TempFull[,1])[1]),1]
SecondAllExceptSpecialisations <- TempFull[c(grep("specialisations",TempFull[,1])[2]:nrow(TempFull)),1]
Tommy <- c(FirstAllExceptSpecialisations,BigSet,SecondAllExceptSpecialisations)
write.table(Tommy,file=file.path(LocationsDir,"TrialP1.xml"),qmethod="escape",quote=F,sep="\t",col.names=F,row.names=F)
cat("Submitting jobs!............")
system(paste("java -jar /lustre/mib-cri/carrol09/MyPipe/workflow-all-1.2-SNAPSHOT.jar --mode=lsf ",file.path(LocationsDir,"TrialP1.xml"),sep=""),wait=TRUE,intern=FALSE)
cat("Jobs Submitted!\n")

sampleSheet <- read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=",",stringsAsFactors=F)
LogFiles <- dir(path="./bamFiles",pattern="*_fileLog.log$")
ProcessedFiles <- dir(path=BamDir,pattern="*_Processed.bam$")
BamsInSampleSheet <-  as.vector(sampleSheet[,"bamFileName"])
for(i in 1:length(BamsInSampleSheet)){
if(any(gsub("_Processed.bam","",ProcessedFiles) %in% gsub(".bam","",BamsInSampleSheet[i]))){
   sampleSheet[i,"Processed_bamFileName"] <- ProcessedFiles[gsub("_Processed.bam","",ProcessedFiles) %in% gsub(".bam","",BamsInSampleSheet[i])]
   LogToRead <- file.path(BamDir,LogFiles[gsub("_fileLog.log","",LogFiles) %in% gsub(".bam","",BamsInSampleSheet[i])])
   DataIn <- read.delim(LogToRead,sep="\t")
   sampleSheet[i,c("Original","delRand","Excluded","Filtered","Unique")] <- DataIn[,c("Mapped","NonRandomChr","IncludedRegions","QC...15","Unique")]
   sampleSheet[i,"DuplicationRate"] <- (as.numeric(sampleSheet[i,"Filtered"])-as.numeric(sampleSheet[i,"Unique"]))/as.numeric(sampleSheet[i,"Filtered"])*100
   sampleSheet[i,"NRF"] <- DataIn[,"Unique"]/DataIn[,"QC...15"]
}else{
	sampleSheet[i,"Processed_bamFileName"] <- "No_Processed_Bam"
	sampleSheet[i,c("Original","delRand","Excluded","Filtered","Unique")] <- rep("No_Information_Available",5)
	sampleSheet[i,"DuplicationRate"] <- "No_Information_Available"
}


}


FragFiles <- dir(path=file.path(WkgDir,"Fragment_Lengths"),pattern="*AllFragLog$")



BamsInSampleSheet <-  as.vector(sampleSheet[,"bamFileName"])
JustBams <- gsub(".bam","",BamsInSampleSheet)
BamsWithFrags <- JustBams[JustBams %in% gsub("_Processed.AllFragLog","",FragFiles)]


for(i in 1:length(BamsInSampleSheet)){
	if(gsub(".bam","",BamsInSampleSheet[i]) %in% BamsWithFrags){
		FragsToRead <- file.path(WkgDir,"Fragment_Lengths",FragFiles[gsub("_Processed.AllFragLog","",FragFiles) %in% gsub(".bam","",BamsInSampleSheet[i])])
		print(FragsToRead)
		DataIn <- as.matrix(read.delim(FragsToRead,sep=" ",comment.char="#"))
		sampleSheet[i,"Sissr_Fragment_Length"] <- DataIn[1,1]		
		sampleSheet[i,"Correlation_Fragment_Length"] <- DataIn[1,2]
		sampleSheet[i,"Coverage_Fragment_Length"] <- DataIn[1,3]
	}else{
		sampleSheet[i,"Sissr_Fragment_Length"] <- "No_Information_Available"		
		sampleSheet[i,"Correlation_Fragment_Length"] <- "No_Information_Available"
		sampleSheet[i,"Coverage_Fragment_Length"] <- "No_Information_Available"	

	}
}

CorrFiles <- dir(path=file.path(WkgDir,"Fragment_Lengths"),pattern="*CorrFragLog$")
BamsWithCorrs <- JustBams[JustBams %in% gsub("_Processed.CorrFragLog","",CorrFiles)]

for(i in 1:length(BamsInSampleSheet)){
	if(gsub(".bam","",BamsInSampleSheet[i]) %in% BamsWithCorrs){
		CorrToRead <- file.path(WkgDir,"Fragment_Lengths",CorrFiles[gsub("_Processed.CorrFragLog","",CorrFiles) %in% gsub(".bam","",BamsInSampleSheet[i])])
		print(CorrToRead)
		DataIn <- as.matrix(read.delim(CorrToRead,sep=" ",comment.char="#"))
		sampleSheet[i,"CC_Fragment_Length"] <- DataIn[1,1]		
		sampleSheet[i,"NSC"] <- DataIn[1,2]
		sampleSheet[i,"RSC"] <- DataIn[1,3]
	}else{
		sampleSheet[i,"CC_Fragment_Length"] <- "No_Information_Available"		
		sampleSheet[i,"NSC"] <- "No_Information_Available"
		sampleSheet[i,"RSC"] <- "No_Information_Available"	
	}
}


write.table(sampleSheet,file=file.path(WkgDir,"SampleSheet.csv"),row.names=F,sep=",")



TempMat <- t(data.matrix(sampleSheet[!sampleSheet[,"Processed_bamFileName"] %in% "No_Processed_Bam",c("Original","delRand","Excluded","Filtered","Unique"),]))
TempTop <- rbind(TempMat[1,]-TempMat[2,],TempMat[2,]-TempMat[3,],TempMat[3,]-TempMat[4,],TempMat[4,]-TempMat[5,],TempMat[5,])[5:1,]



png(file.path(WkgDir,"bamFiles","ReadsPlot.png"))
barplot(TempTop,names=sampleSheet[!sampleSheet[i,"Processed_bamFileName"] %in% "No_Processed_Bam","SampleID"],horiz=T,las=1,legend.text=c("Pass Filter","Duplicated","MapQ < 15","Excluded Regions","Random_Contigs"))
dev.off()




write.table("Complete",file.path(LocationsDir,paste(Args[1],"_MainBamProcess.txt",sep="")),col.names=T,row.names=F,sep=",",quote=F)


