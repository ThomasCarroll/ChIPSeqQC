getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))

library(raster)
Args <- commandArgs(trailingOnly = TRUE)
ConfigDirectory <- Args[2]
ConfigFile <- readIniFile(file.path(ConfigDirectory,"config.ini"))
TempDir <- ConfigFile[ConfigFile[,2] %in% "tempdirectory",3]
TempFull <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/ReportMeta.xml",header=F,stringsAsFactors=F,quote="")


TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1] <- gsub("WorkingDirectory_ToChange",ConfigFile[ConfigFile[,2] %in% "workingdirectory",3],TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1] <- gsub("TempDirectory_ToChange",ConfigFile[ConfigFile[,2] %in% "tempdirectory",3],TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1])

write.table(TempFull,file=file.path(TempDir,"ChIPReportMeta.xml"),qmethod="escape",quote=F,sep="\t",col.names=F,row.names=F)
cat("Submitting jobs!............")
system(paste("java -jar /lustre/mib-cri/carrol09/MyPipe/workflow-all-1.2-SNAPSHOT.jar --mode=lsf ",file.path(TempDir,"ChIPReportMeta.xml"),sep=""),wait=FALSE,intern=FALSE)
cat("Jobs Submitted!\n")
write.table("Complete",file.path(ConfigDirectory,paste(Args[1],"_MainReportProcess.txt",sep="")),col.names=T,row.names=F,sep=",",quote=F)