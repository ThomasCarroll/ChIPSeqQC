Arguments <- commandArgs(trailingOnly = T)
getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
fileIn <- Arguments[1]
DirectoryOut <- Arguments[2]

library(GenomicRanges)
library(htSeqTools)
library(Rsamtools)

#file="/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/AnotherVerison/20111109_RossAdams_DN_HNF1bChIP/bamFiles/SLX-4499.739.s_4.bwa.homo_sapiens_Processed.bam"
#DirectoryOut <- "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/AnotherVerison/20111109_RossAdams_DN_HNF1bChIP/"

Temp <- scanBamHeader(fileIn)
AllSeqnames <- Temp[[1]]$targets
LongestChr <- AllSeqnames[which.max(AllSeqnames)]
LCName <- names(LongestChr)
LCLength <- unique(LongestChr)

which <- GRanges(seqnames=LCName,ranges=IRanges(1,LCLength))
param <- ScanBamParam(which=which)


Con_input<-readGappedAlignments(fileIn,param=param)
Con_input <- granges(Con_input)

#samples<-RangedDataList(H3k4_1=as(Con_input,"RangedData"),Pol2=as(Second_input,"RangedData"),Input=as(Third_input,"RangedData"),Myc=as(Fourth_input,"RangedData"))
samples<-RangedDataList(Example=as(Con_input,"RangedData"))
#cmds1 <- cmds(samples,k=2)
ssdOfSample <- ssdCoverage(samples)
giniOfSample <- giniCoverage(samples)
save(ssdOfSample,giniOfSample,file=paste(DirectoryOut,"/",gsub("/.*/","",gsub(".bam","",fileIn)),"_CoverageDispersion.RData",sep=""))
#png(paste(DirectoryOut,"/",gsub("/.*/","",gsub(".bam","",fileIn)),"_GiniCoverage.png",sep=""))
#giniCoverage(samples[[1]],mc.cores=1,mk.plot=TRUE)
#dev.off()



