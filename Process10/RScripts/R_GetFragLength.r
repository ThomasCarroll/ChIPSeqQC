Arguments <- commandArgs(trailingOnly = T)
getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
file <- Arguments[1]
DirectoryOut <- Arguments[2]

library(chipseq)
library(GenomicRanges)

NewName <- gsub(".bam","",gsub("/.*/","",file))


Temp <- scanBamHeader(file)
AllSeqnames <- Temp[[1]]$targets
LongestChr <- AllSeqnames[which.max(AllSeqnames)]
LCName <- names(LongestChr)
LCLength <- unique(LongestChr)

which <- GRanges(seqnames=LCName,ranges=IRanges(1,LCLength))
param <- ScanBamParam(which=which)


BamBnd <-  readGappedAlignments(file,param=param)
if(length(BamBnd) > 1000000){
GrangesAlign <- granges(BamBnd)



GrangesAlign2 <- GrangesAlign

ReadLength <- round(median(width(GrangesAlign2)))

 
Temp2 <- list(start(GrangesAlign2[strand(GrangesAlign2) == "+"]),start(GrangesAlign2[strand(GrangesAlign2) == "-"]))
names(Temp2) <- c("+","-")



fraglenSissr <- median(estimate.mean.fraglen(GrangesAlign2))
#fraglenCorr <- median(estimate.mean.fraglen(GrangesAlign,method="correlation"))
fraglenCorr <- NA
fraglenCov <- median(estimate.mean.fraglen(GrangesAlign2,method="coverage"))

CovPlot <- basesCovered(Temp2, shift = seq(0, 300, 5),seqLen=36)

#png(paste(file.path(DirectoryOut,NewName),"Frag_Cov.png",sep="_"))
#plot(CovPlot,type="l",col="red")
#dev.off()



AlFragLens <- cbind(fraglenSissr,fraglenCorr,fraglenCov)
colnames(AlFragLens) <- c("Sissr","Correlation","Coverage")
write.table(AlFragLens,paste(file.path(DirectoryOut,NewName),".AllFragLog",sep=""),row.names=FALSE)
write.table(CovPlot,paste(file.path(DirectoryOut,NewName),".FragCovLog",sep=""),row.names=FALSE,col.names=F)
}else{
AlFragLens <- matrix(c("Too_few_Reads_To_Calculate","Too_few_Reads_To_Calculate","Too_few_Reads_To_Calculate"),ncol=3,nrow=1)
colnames(AlFragLens) <- c("Sissr","Correlation","Coverage")
write.table(AlFragLens,paste(file.path(DirectoryOut,NewName),".AllFragLog",sep=""),row.names=FALSE)
}
