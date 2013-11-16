getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
source("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/RScripts/Workflow_Functions.r")

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

if(all(c("Macs_name","Sicer_Name","TPICS_Name") %in% colnames(sampleSheet))){

SampleNames <- sampleSheet[!is.na(sampleSheet[,"Macs_name"]) & !is.na(sampleSheet[,"Sicer_Name"]) & !is.na(sampleSheet[,"TPICS_Name"]),"SampleName"]







TempFull <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/AcrossCallersMeta.xml",header=F,stringsAsFactors=F,quote="")
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
  AllSets[[i]][grep("samplesheet_ToChange",AllSets[[i]])] <- gsub("samplesheet_ToChange",file.path(getwd(),"SampleSheet.csv"),AllSets[[i]][grep("samplesheet_ToChange",AllSets[[i]])])

  BigSet <- c(BigSet,AllSets[[i]])
}

FirstAllExceptSpecialisations <- TempFull[c(1:grep("specialisations",TempFull[,1])[1]),1]
SecondAllExceptSpecialisations <- TempFull[c(grep("specialisations",TempFull[,1])[2]:nrow(TempFull)),1]
Tommy <- c(FirstAllExceptSpecialisations,BigSet,SecondAllExceptSpecialisations)
write.table(Tommy,file=file.path(LocationsDir,"RunAcrossCallers.xml"),qmethod="escape",quote=F,sep="\t",col.names=F,row.names=F)
cat("Submitting jobs!............")
system(paste("java -Xmx1G -jar /lustre/mib-cri/carrol09/MyPipe/workflow-all-1.2-SNAPSHOT.jar --mode=lsf ",file.path(LocationsDir,"RunAcrossCallers.xml"),sep=""),wait=TRUE,intern=FALSE)



#sampleSheet <- read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=",",stringsAsFactors=F)
sampleSheet <- ReadAndLock(file.path(WkgDir,"SampleSheet.csv"),WkdDir,SAF=F,napTime=5)
for(i in 1:nrow(sampleSheet)){
	SampleToLookFor <- sampleSheet[i,"SampleName"]
	PeakAcrossFile <- dir(path=file.path(WkgDir,"Peaks","PeakProfiles"),pattern=paste(SampleToLookFor,"_AcrossPeakCallers.txt",sep=""),full.names=T)
	if(length(PeakAcrossFile) > 0){
		Temp <- read.delim(PeakAcrossFile,sep="\t",h=T)
		sampleSheet[i,"Macs_In_TPICs"] <- as.vector(Temp[1,2])
		sampleSheet[i,"TPICs_In_Macs"] <- as.vector(Temp[2,2])
		sampleSheet[i,"Macs_In_Sicer"]  <- as.vector(Temp[3,2])
		sampleSheet[i,"Sicer_In_Macs"]  <- as.vector(Temp[4,2])
		sampleSheet[i,"TPICs_In_Sicer"] <- as.vector(Temp[5,2])
		sampleSheet[i,"Sicer_In_TPICs"] <- as.vector(Temp[6,2])
		sampleSheet[i,"JacardIndex_TPICS_And_Macs"] <- as.vector(Temp[7,2])
		sampleSheet[i,"JacardIndex_Sicer_And_Macs"] <- as.vector(Temp[8,2])
		sampleSheet[i,"JacardIndex_TPICs_And_Sicer"] <- as.vector(Temp[9,2])
 	}else{
		sampleSheet[i,"Macs_In_TPICs"] <- ""
		sampleSheet[i,"TPICs_In_Macs"] <- ""
		sampleSheet[i,"Macs_In_Sicer"]  <- ""
		sampleSheet[i,"Sicer_In_Macs"]  <- ""
		sampleSheet[i,"TPICs_In_Sicer"] <- ""
		sampleSheet[i,"Sicer_In_TPICs"] <- ""
		sampleSheet[i,"JacardIndex_TPICS_And_Macs"] <- ""
		sampleSheet[i,"JacardIndex_Sicer_And_Macs"] <- ""
		sampleSheet[i,"JacardIndex_TPICs_And_Sicer"] <- ""
	}
}


#write.table(sampleSheet,file=file.path(WkgDir,"SampleSheet.csv"),row.names=F,sep=",")
WriteAndUnlock(sampleSheet,file.path(WkgDir,"SampleSheet.csv"))	

}
write.table("Complete",file.path(LocationsDir,paste(Args[1],"_MainAcrossPeaksProcess.txt",sep="")),col.names=T,row.names=F,sep=",",quote=F)
