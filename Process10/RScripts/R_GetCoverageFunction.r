Arguments <- commandArgs(trailingOnly = T)
getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
file <- Arguments[1]
DirectoryOut <- Arguments[2]

NewName <- gsub(".bam","",gsub("/.*/","",file))
BamBnd <-  readGappedAlignments(file)
GrangesAlign <- granges(BamBnd)
#GrangesAlign <- GrangesAlign[seqnames(GrangesAlign) %in% c("chr1")]
ReadLength <- round(median(width(GrangesAlign)))


fraglenSissr <- median(estimate.mean.fraglen(GrangesAlign))
fraglenCorr <- median(estimate.mean.fraglen(GrangesAlign,method="correlation"))
fraglenCov <- median(estimate.mean.fraglen(GrangesAlign,method="coverage"))

CovPlot <- basesCovered(Temp2, shift = seq(0, 300, 5),seqLen=36)
CorrPlot <- densityCorr(Temp2, shift = seq(0, 300, 5), center = FALSE,seqLen=36)

png(paste(file.path(DirectoryOut,NewName),"Frag_Cov.png",sep="_"))
plot(CovPlot,type="l",col="red")
dev.off()

png(paste(file.path(DirectoryOut,NewName),"Frag_Corr.png",sep="_"))
plot(CorrPlot,type="l",col="red")
dev.off()


AlFragLens <- cbind(fraglenSissr,fraglenCorr,fraglenCov)
colnames(AlFragLens) <- c("Sissr","Correlation","Coverage")
write.table(,paste(file.path(DirectoryOut,NewName),".AllFragLog",sep=""))
write.table(CovPlot,paste(file.path(DirectoryOut,NewName),".FragCovLog",sep=""),row.names=FALSE,col.names=F)
write.table(CorrPlot,paste(file.path(DirectoryOut,NewName),"FragCorrLog",sep=""),row.names=FALSE,col.names=F)
