Arguments <- commandArgs(trailingOnly = T)
getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))

#WkgDir <- Arguments[1]
WkgDir <- getwd()
sampleSheet <- read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=",",stringsAsFactors=F)
LogFiles <- dir(path="./bamFiles",pattern="*_fileLog.log$")
ProcessedFiles <- dir(path="./bamFiles",pattern="*_Processed.bam$")
BamsInSampleSheet <-  as.vector(sampleSheet[,"bamFileName"])
for(i in 1:length(BamsInSampleSheet)){
if(any(gsub("_Processed.bam","",ProcessedFiles) %in% gsub(".bam","",BamsInSampleSheet[i]))){
   sampleSheet[i,"Processed_bamFileName"] <- ProcessedFiles[gsub("_Processed.bam","",ProcessedFiles) %in% gsub(".bam","",BamsInSampleSheet[i])]
   LogToRead <- file.path(WkgDir,"bamFiles",LogFiles[gsub("_fileLog.log","",LogFiles) %in% gsub(".bam","",BamsInSampleSheet[i])])
   DataIn <- read.delim(LogToRead,sep="\t")
   sampleSheet[i,c("Original","delRand","Excluded","Filtered","Unique")] <- DataIn[,c("Mapped","NonRandomChr","IncludedRegions","QC...15","Unique")]
   sampleSheet[i,"DuplicationRate"] <- (as.numeric(sampleSheet[i,"Filtered"])-as.numeric(sampleSheet[i,"Unique"]))/as.numeric(sampleSheet[i,"Filtered"])
}else{
	sampleSheet[i,"Processed_bamFileName"] <- "No_Processed_Bam"
	sampleSheet[i,c("Original","delRand","Excluded","Filtered","Unique")] <- rep("No_Information_Available",5)
	sampleSheet[i,"DuplicationRate"] <- "No_Information_Available"
}


}
write.table(sampleSheet,file=file.path(WkgDir,"SampleSheet.csv"),row.names=F,sep=",")



TempMat <- t(data.matrix(sampleSheet[!sampleSheet[,"Processed_bamFileName"] %in% "No_Processed_Bam",c("Original","delRand","Excluded","Filtered","Unique"),]))
TempTop <- rbind(TempMat[1,]-TempMat[2,],TempMat[2,]-TempMat[3,],TempMat[3,]-TempMat[4,],TempMat[4,]-TempMat[5,],TempMat[5,])[5:1,]



png(file.path(WkgDir,"bamFiles","ReadsPlot.png"))
barplot(TempTop,names=sampleSheet[!sampleSheet[i,"Processed_bamFileName"] %in% "No_Processed_Bam","SampleID"],horiz=T,las=1,legend.text=c("Pass Filter","Duplicated","MapQ < 15","Excluded Regions","Random_Contigs"))
dev.off()

