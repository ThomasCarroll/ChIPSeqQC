Arguments <- commandArgs(trailingOnly = T)
getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
WkgDir <- Arguments[1]

#WkgDir <- "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet"

sampleSheet <- as.matrix(read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=","))


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

write.table(sampleSheet,file=file.path(WkgDir,"SampleSheet.csv"),row.names=F,sep=",")
