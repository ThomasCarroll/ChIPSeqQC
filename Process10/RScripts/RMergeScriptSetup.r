getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))

Args <- commandArgs(trailingOnly = TRUE)
library(raster)
ConfigFile <- readIniFile(file.path(Args[2],"config.ini"))
ConfigFile <- readIniFile(file.path(getwd(),"Temp","config.ini"))
LocationsDir <- ConfigFile[ConfigFile[,2] %in% "tempdirectory",3]
WkgDir<- ConfigFile[ConfigFile[,2] %in% "workingdirectory",3]
BamDir<- ConfigFile[ConfigFile[,2] %in% "bamdirectory",3]
FQDir<- ConfigFile[ConfigFile[,2] %in% "fastqdirectory",3]


sampleSheet <- read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=",")
UniquedMergedFiles <-  as.vector(unique(na.omit(sampleSheet[,"toMerge"])))
AllSamplesPresent <- as.vector(unique(na.omit(sampleSheet[,"GenomicsID"])))
UniquedMergedFiles <- UniquedMergedFiles[!UniquedMergedFiles %in% AllSamplesPresent]
if(length(UniquedMergedFiles) != 0){
  InputLists <- vector("character",length=length(UniquedMergedFiles))
  OutputNames <- UniquedMergedFiles
  for(i in 1:length(UniquedMergedFiles)){
      ToBeMerged <-  as.vector(sampleSheet[sampleSheet[,"toMerge"] %in% UniquedMergedFiles[i],"bamFileName"])
	ToBeMerged <- ToBeMerged[!is.na(ToBeMerged)]
      InputLists[i] <- paste(paste("INPUT=",file.path(BamDir,ToBeMerged),sep=""),collapse=" ")
  }
UniquedMergedFiles <- UniquedMergedFiles[!gsub("INPUT=","",InputLists) == ""]
TempFull <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/MergeMeta.xml",header=F,stringsAsFactors=F,quote="")
TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1] <- gsub("WorkingDirectory_ToChange",WkgDir,TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1])
SpecialisationSet <- TempFull[grep("specialisation identifier",TempFull[,1]):max(grep("specialisation>",TempFull[,1])),1]
AllSets <- rep(list(SpecialisationSet),length(UniquedMergedFiles))
BigSet <- vector("character")
for(i in 1:length(AllSets)){
  AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])] <- gsub("Indentifier_ToChange",getRandString(len=6),AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("Input_To_Fill",AllSets[[i]])] <- gsub("Input_To_Fill",InputLists[i],AllSets[[i]][grep("Input_To_Fill",AllSets[[i]])])
  AllSets[[i]][grep("Output_To_Fill",AllSets[[i]])] <- gsub("Output_To_Fill",OutputNames[i],AllSets[[i]][grep("Output_To_Fill",AllSets[[i]])])
  BigSet <- c(BigSet,AllSets[[i]])
}

FirstAllExceptSpecialisations <- TempFull[c(1:grep("specialisations",TempFull[,1])[1]),1]
SecondAllExceptSpecialisations <- TempFull[c(grep("specialisations",TempFull[,1])[2]:nrow(TempFull)),1]
Tommy <- c(FirstAllExceptSpecialisations,BigSet,SecondAllExceptSpecialisations)
write.table(Tommy,file=file.path(WkgDir,"TrialMerge.xml"),qmethod="escape",quote=F,sep="\t",col.names=F,row.names=F)
cat("Submitting jobs!............")
system(paste("java -jar /lustre/mib-cri/carrol09/MyPipe/workflow-all-1.2-SNAPSHOT.jar --mode=lsf ",file.path(WkgDir,"TrialMerge.xml"),sep=""),wait=TRUE,intern=FALSE)
cat("Jobs Submitted!\n")
NewSampleSheet <- matrix(ncol=ncol(sampleSheet),nrow=length(OutputNames))
colnames(NewSampleSheet) <- colnames(sampleSheet)
NewSampleSheet <- as.matrix(rbind(sampleSheet,NewSampleSheet))
	for(i in 1:length(OutputNames)){
		NewSampleSheet[nrow(sampleSheet)+i,"GenomicsID"] <- OutputNames[i]
		NewSampleSheet[nrow(sampleSheet)+i,"SampleName"] <- OutputNames[i]
		NewSampleSheet[nrow(sampleSheet)+i,"bamFileName"] <- paste(OutputNames[i],".bwa.bam",sep="")
	}
write.table(NewSampleSheet,"SampleSheet.csv",sep=",",row.names=F,quote=F)
}
write.table("Complete",file.path(LocationsDir,paste(Args[1],"_MainMerge.txt",sep="")),col.names=T,row.names=F,sep=",",quote=F)


