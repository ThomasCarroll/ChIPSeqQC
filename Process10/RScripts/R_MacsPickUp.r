Arguments <- commandArgs(trailingOnly = T)
getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
WkgDir <- Arguments[1]

#WkgDir <- "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet"

sampleSheet <- read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=",")


MacsFiles <- dir(path=file.path(WkgDir,"Peaks","Macs_Peaks"),pattern="*peaks.xls$")
PositivePeaks <- MacsFiles[-grep("negative",MacsFiles)]


BamsInSampleSheet <-  as.vector(sampleSheet[,"bamFileName"])
JustBams <- gsub(".bam","",BamsInSampleSheet)
BamsWithPeaks <- JustBams[JustBams %in% gsub("_Processed_peaks.xls","",MacsFiles)]


for(i in 1:length(BamsInSampleSheet)){
	if(gsub(".bam","",BamsInSampleSheet[i]) %in% BamsWithPeaks){
		MacsToRead <- file.path(WkgDir,"Peaks","Macs_Peaks",MacsFiles[gsub("_Processed_peaks.xls","",MacsFiles) %in% gsub(".bam","",BamsInSampleSheet[i])])
		print(MacsToRead)
		DataIn <- read.delim(MacsToRead,sep="\t",comment.char="#")
		sampleSheet[i,"macsPeaks"] <- nrow(DataIn)		
		sampleSheet[i,"macs_name"] <- MacsToRead
		print(nrow(DataIn))
	}
}

write.table(sampleSheet,file=file.path(WkgDir,"SampleSheet.csv"),row.names=F,sep=",")
