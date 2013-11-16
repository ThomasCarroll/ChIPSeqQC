Arguments <- commandArgs(trailingOnly = T)
getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
file <- Arguments[1]
DirectoryOut <- Arguments[2]
#ReadLength <- Arguments[3]
library(chipseq)
library(GenomicRanges)
library(snow)
library(spp)

BamBnd <-  readGappedAlignments(file)
NewName <- gsub(".bam","",gsub("/.*/","",file))
if(length(BamBnd) > 1000000){

GrangesAlign <- granges(BamBnd)
ReadLength <- round(median(width(GrangesAlign)))


tempo <-  read.bam.tags(file)
binding.characteristics <- get.binding.characteristics(tempo,srange=c(0,400),bin=2,accept.all.tags=T,remove.tag.anomalies = F)
CorMaxPos <- binding.characteristics$peak$x
CorMin <- min(binding.characteristics$cross.correlation[,2])
CorMax <- max(binding.characteristics$cross.correlation[,2])
CorReadLength <- binding.characteristics$cross.correlation[which.min(abs(binding.characteristics$cross.correlation[,1]-ReadLength)),2]
write.table(binding.characteristics$cross.correlation,paste(file.path(DirectoryOut,NewName),".CorrFragfull",sep=""),row.names=FALSE)
NSC <- CorMax/CorMin
ReadLenghSC <- CorReadLength/CorMin
RSC <- NSC/ReadLenghSC
  ToPrint <- cbind(CorMaxPos,NSC,RSC)
}else{
  ToPrint <- cbind("Too_few_Reads_To_Calculate","Too_few_Reads_To_Calculate","Too_few_Reads_To_Calculate")
}
write.table(ToPrint,paste(file.path(DirectoryOut,NewName),".CorrFragLog",sep=""),row.names=FALSE)

#GetCorrForChr <- function(Chromosome,GrangesAlign){
#  GrangesAlign2 <- GrangesAlign[seqnames(GrangesAlign) %in% Chromosome]
#  Temp2 <- list(start(GrangesAlign2[strand(GrangesAlign2) == "+"]),start(GrangesAlign2[strand(GrangesAlign2) == "-"]))
#  names(Temp2) <- c("+","-")
#  Tommy <- densityCorr(Temp2, shift = seq(0, 400, 5), center = FALSE,width = ReadLength *2L, seqLen=ReadLength)
#  rm(Temp2)
#  return(Tommy)
#}


#NewName <- gsub(".bam","",gsub("/.*/","",file))
#BamBnd <-  readGappedAlignments(file)

#Chrs <- list("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")


#if(length(BamBnd) > 1000000){
 # GrangesAlign <- granges(BamBnd)
 # ReadLength <- round(median(width(GrangesAlign)))
  #GrangesAlign <- GrangesAlign[seqnames(GrangesAlign) %in% unlist(Chrs)]
 # ResCor <- mclapply(Chrs,GetCorrForChr,GrangesAlign=GrangesAlign)
#  GetCorrForChr(GrangesAlign,Chrs[[1]])
#}

#if(length(BamBnd) > 1000000){
#GrangesAlign <- granges(BamBnd)
#ReadLength <- round(median(width(GrangesAlign)))
#Tommy <- vector("list",length=length(unique(seqnames(GrangesAlign))))
#for(i in 1:length(unique(seqnames(GrangesAlign)))){
#  GrangesAlign2 <- GrangesAlign[seqnames(GrangesAlign) %in% unique(seqnames(GrangesAlign))[[i]]]
#  Temp2 <- list(start(GrangesAlign2[strand(GrangesAlign2) == "+"]),start(GrangesAlign2[strand(GrangesAlign2) == "-"]))
#  names(Temp2) <- c("+","-")
#  Tommy[[i]] <- densityCorr(Temp2, shift = seq(0, 300, 5), center = FALSE,width = ReadLength *2L, seqLen=ReadLength)
#  print(i)
#  rm(Temp2)
#}
#LongCorVec <- rowMeans(Tommy)


#AlFragLens <- cbind(fraglenCorr)
#colnames(AlFragLens) <- c("Correlation")
#write.table(AlFragLens,paste(file.path(DirectoryOut,"Fragment_Lengths",NewName),".AllFragLog",sep=""),row.names=FALSE)
#write.table(CovPlot,paste(file.path(DirectoryOut,"Fragment_Lengths",NewName),".FragCovLog",sep=""),row.names=FALSE,col.names=F)
#}else{
#AlFragLens <- matrix(c("Too_few_Reads_To_Calculate","Too_few_Reads_To_Calculate","Too_few_Reads_To_Calculate"),ncol=3,nrow=1)
#colnames(AlFragLens) <- c("Sissr","Correlation","Coverage")
#write.table(AlFragLens,paste(file.path(DirectoryOut,"Fragment_Lengths",NewName),".AllFragLog",sep=""),row.names=FALSE)
#}