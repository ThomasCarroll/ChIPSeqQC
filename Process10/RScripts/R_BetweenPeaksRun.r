getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
Args <- commandArgs(trailingOnly = TRUE)
library(raster)

source("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/RScripts/Workflow_Functions.r")

#ConfigFile <- readIniFile(file.path(Args[2],"config.ini"))
ConfigFile <- readIniFile(file.path(getwd(),"Temp","config.ini"))
LocationsDir <- ConfigFile[ConfigFile[,2] %in% "tempdirectory",3]
WkgDir<- ConfigFile[ConfigFile[,2] %in% "workingdirectory",3]
BamDir<- ConfigFile[ConfigFile[,2] %in% "bamdirectory",3]
Genome <- ConfigFile[ConfigFile[,2] %in% "genome",3]
GenomeFileOptions <-  ConfigFile[ConfigFile[,1] %in% "Genomes",]
GenomeFile <- GenomeFileOptions[GenomeFileOptions[,2] %in% tolower(Genome),3]

sampleSheet <- read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=",")
if(all(c("Macs_name") %in% colnames(sampleSheet))){

Peaks <- as.vector(sampleSheet[!is.na(sampleSheet[,"Macs_name"]),"Macs_name"])
samples <- as.vector(sampleSheet[!is.na(sampleSheet[,"Macs_name"]),"SampleName"])
SampleNames <- vector("character",length=length(Peaks)*length(Peaks))
SampleNames2 <- vector("character",length=length(Peaks)*length(Peaks))
OutNames <- vector("character",length=length(Peaks)*length(Peaks))
OutNames2 <- vector("character",length=length(Peaks)*length(Peaks))

J <- 1
for(i in 1:length(Peaks)){
  for(k in 1:length(Peaks)){
    SampleNames[J] <- Peaks[i]
    SampleNames2[J] <- Peaks[k]
    OutNames[J] <- samples[i]
    OutNames2[J] <- samples[k]
    J <- J+1
  }
}

OutNameFull <- paste(OutNames,"_Vs_",OutNames2,sep="")

TempFull <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/BetweenPeakSetsMeta.xml",header=F,stringsAsFactors=F,quote="")
TempFull[grep("RunDirectory_ToChange",TempFull[,1]),1] <- gsub("RunDirectory_ToChange",getwd(),TempFull[grep("RunDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("BamDirectory_ToChange",TempFull[,1]),1] <- gsub("BamDirectory_ToChange",BamDir,TempFull[grep("BamDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1] <- gsub("WorkingDirectory_ToChange",WkgDir,TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1] <- gsub("TempDirectory_ToChange",LocationsDir,TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1])


SpecialisationSet <- TempFull[grep("specialisation identifier",TempFull[,1]):max(grep("specialisation>",TempFull[,1])),1]
AllSets <- rep(list(SpecialisationSet),length(SampleNames ))
BigSet <- vector("character")
for(i in 1:length(AllSets)){
  AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])] <- gsub("Indentifier_ToChange",getRandString(len=6),AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("Test_ToChange",AllSets[[i]])] <- gsub("Test_ToChange",SampleNames[i],AllSets[[i]][grep("Test_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("Test2_ToChange",AllSets[[i]])] <- gsub("Test2_ToChange",SampleNames2[i],AllSets[[i]][grep("Test2_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("OutName_ToChange",AllSets[[i]])] <- gsub("OutName_ToChange",OutNameFull[i],AllSets[[i]][grep("OutName_ToChange",AllSets[[i]])])
  BigSet <- c(BigSet,AllSets[[i]])
}

FirstAllExceptSpecialisations <- TempFull[c(1:grep("specialisations",TempFull[,1])[1]),1]
SecondAllExceptSpecialisations <- TempFull[c(grep("specialisations",TempFull[,1])[2]:nrow(TempFull)),1]
Tommy <- c(FirstAllExceptSpecialisations,BigSet,SecondAllExceptSpecialisations)
write.table(Tommy,file=file.path(LocationsDir,"RunBetweenMacs.xml"),qmethod="escape",quote=F,sep="\t",col.names=F,row.names=F)
cat("Submitting jobs!............")
system(paste("java -Xmx1G -jar /lustre/mib-cri/carrol09/MyPipe/workflow-all-1.2-SNAPSHOT.jar --mode=lsf ",file.path(LocationsDir,"RunBetweenMacs.xml"),sep=""),wait=TRUE,intern=FALSE)
}
write.table("Complete",file.path(LocationsDir,paste(Args[1],"_MainBetweenPeaksProcess.txt",sep="")),col.names=T,row.names=F,sep=",",quote=F)


