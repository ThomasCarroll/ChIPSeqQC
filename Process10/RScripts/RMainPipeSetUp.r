getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))

library(raster)
Args <- commandArgs(trailingOnly = TRUE)
ConfigDirectory <- Args[1]
ConfigFile <- readIniFile(file.path(ConfigDirectory,"config.ini"))
TempDir <- ConfigFile[ConfigFile[,2] %in% "tempdirectory",3]
TempFull <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/MainPipe.xml",header=F,stringsAsFactors=F,quote="")


TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1] <- gsub("WorkingDirectory_ToChange",ConfigFile[ConfigFile[,2] %in% "workingdirectory",3],TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1] <- gsub("TempDirectory_ToChange",ConfigFile[ConfigFile[,2] %in% "tempdirectory",3],TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("FQDirectory_ToChange",TempFull[,1]),1] <- gsub("FQDirectory_ToChange",ConfigFile[ConfigFile[,2] %in% "fastqdirectory",3],TempFull[grep("FQDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("BamDirectory_ToChange",TempFull[,1]),1] <- gsub("BamDirectory_ToChange",ConfigFile[ConfigFile[,2] %in% "bamdirectory",3],TempFull[grep("BamDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("ConFig_ToChange",TempFull[,1]),1] <- gsub("ConFig_ToChange",file.path(ConfigDirectory,"config.ini"),TempFull[grep("ConFig_ToChange",TempFull[,1]),1])
TempFull[grep("PathwayTrackerID_To_Fill",TempFull[,1]),1] <- gsub("PathwayTrackerID_To_Fill",paste("PathID",getRandString(),sep=""),TempFull[grep("PathwayTrackerID_To_Fill",TempFull[,1]),1])
TempFull[grep("Indentifier_ToChange",TempFull[,1]),1] <- gsub("Indentifier_ToChange",paste("SpecialID",getRandString(),sep=""),TempFull[grep("Indentifier_ToChange",TempFull[,1]),1])
TempFull[grep("TemporaryDirectory_ToChange",TempFull[,1]),1] <- gsub("TemporaryDirectory_ToChange",ConfigFile[ConfigFile[,2] %in% "tempdirectory",3],TempFull[grep("TemporaryDirectory_ToChange",TempFull[,1]),1])

write.table(TempFull,file=file.path(TempDir,"MainPipe.xml"),qmethod="escape",quote=F,sep="\t",col.names=F,row.names=F)
cat("Submitting jobs!............")
system(paste("java -jar /lustre/mib-cri/carrol09/MyPipe/workflow-all-1.2-SNAPSHOT.jar --mode=lsf ",file.path(TempDir,"MainPipe.xml"),sep=""),wait=FALSE,intern=FALSE)
cat("Jobs Submitted!\n")
