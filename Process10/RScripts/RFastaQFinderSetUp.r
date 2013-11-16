getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))


Args <- commandArgs(trailingOnly = TRUE)
library(raster)
ConfigFile <- readIniFile(file.path(Args[2],"config.ini"))
#ConfigFile <- readIniFile(file.path(getwd(),"Temp","config.ini"))
LocationsDir <- ConfigFile[ConfigFile[,2] %in% "tempdirectory",3]
WkgDir<- ConfigFile[ConfigFile[,2] %in% "workingdirectory",3]
BamDir<- ConfigFile[ConfigFile[,2] %in% "bamdirectory",3]
 
TempFull <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/FQLookUpMeta.xml",header=F,stringsAsFactors=F,quote="")
TempFull[grep("RunDirectory_ToChange",TempFull[,1]),1] <- gsub("RunDirectory_ToChange",getwd(),TempFull[grep("RunDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("BamDirectory_ToChange",TempFull[,1]),1] <- gsub("BamDirectory_ToChange",BamDir,TempFull[grep("BamDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1] <- gsub("WorkingDirectory_ToChange",WkgDir,TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1] <- gsub("TempDirectory_ToChange",LocationsDir,TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("TemporaryDirectory_ToChange",TempFull[,1]),1] <- gsub("TemporaryDirectory_ToChange",LocationsDir,TempFull[grep("TemporaryDirectory_ToChange",TempFull[,1]),1])

if(file.exists(file.path(LocationsDir,"ActualFQLocations.txt"))){
unlink(file.path(LocationsDir,"ActualFQLocations.txt"))
}

write.table(TempFull,file=file.path(LocationsDir,"FQLookUp.xml"),qmethod="escape",quote=F,sep="\t",col.names=F,row.names=F)
cat("Submitting jobs!............")
system(paste("java -jar /lustre/mib-cri/carrol09/MyPipe/workflow-all-1.2-SNAPSHOT.jar --mode=lsf ",file.path(LocationsDir,"FQLookUp.xml"),sep=""),wait=TRUE,intern=FALSE)
cat("Jobs Submitted!\n")
write.table("Complete",file.path(LocationsDir,paste(Args[1],"_MainFQLookProcess.txt",sep="")),col.names=T,row.names=F,sep=",",quote=F)

