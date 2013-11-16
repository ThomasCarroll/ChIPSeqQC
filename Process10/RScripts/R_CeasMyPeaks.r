#PeakAnnoFile <- "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/20110610_MohammadH_JC_Greb1ChIP/Peaks/Macs_Peaks/SLX-2474.412.s_3.bwa.RealignedGRCh37_Processed_peaks.xls"
Args <- commandArgs(trailingOnly = TRUE)
PeakAnnoFile <- Args[1]
OutFileCounts <- Args[2]
OutFileProfile <- Args[3]
PeakAnnoFile <- gsub(".bed","_Annotated.xls",PeakAnnoFile)
PeakAnno <- read.delim(PeakAnnoFile,sep="\t")

IntergenicPeaks <- PeakAnno[PeakAnno[,"Feature"] == "Intergenic",]

PositiveGenes <- PeakAnno[PeakAnno[,"strand"] == "+",]
NegativeGenes <- PeakAnno[PeakAnno[,"strand"] == "-",]
DistanceToPart1 <- (as.vector(PositiveGenes[,"Distance.to.5..end.of.Feature"])+2000)
DistanceToPart2 <- (as.vector(NegativeGenes[,"Distance.to.3..end.of.Feature"])-2000)*-1

Temp <- cbind(rbind(PositiveGenes,NegativeGenes),c(DistanceToPart1,DistanceToPart2))
Temp <- Temp[order(abs(Temp[,"c(DistanceToPart1, DistanceToPart2)"]),decreasing=F),]
TempTrimmed <- Temp[match(unique(Temp[,"Peak_V4"]),Temp[,"Peak_V4"]),]
 NewB <- seq(-10000,2000,by=50)

Tom <- hist(TempTrimmed[TempTrimmed[,ncol(TempTrimmed)] > -10000 & TempTrimmed[,ncol(TempTrimmed)] < 2000,ncol(TempTrimmed)],breaks=NewB)
MainPeaks <- TempTrimmed[TempTrimmed[,ncol(TempTrimmed)] > -10000 & TempTrimmed[,ncol(TempTrimmed)] < 2000,"Peak_V4"]
Counts <- Tom$counts
Breaks <- Tom$breaks
midpoints <- Tom$mids
CountsInSections <- matrix(nrow=7,ncol=3)
CountsInSections[1:7,1] <- c("10000upstream_To_5000upstream","5000upstream_To_2000upstream","2000upstream_To_500upstream","500upstream_To_500downstream","500downastream_To_2000downstream","2000downstream_To_Gene_End","Intergenic")
CountsInSections[1:5,2] <- c(sum(Counts[1:100]),sum(Counts[101:160]),sum(Counts[161:190]),sum(Counts[191:210]),sum(Counts[211:240]))
FurtherCount <- TempTrimmed[!(TempTrimmed[,ncol(TempTrimmed)] > -10000 & TempTrimmed[,ncol(TempTrimmed)] < 2000),]
#OverlappingGene <- FurtherCount[FurtherCount[,"Feature"] == "2000Extended","Peak_V4"]
#Intergenic<- FurtherCount[FurtherCount[,"Feature"] == "Intergenic","Peak_V4"]

CountsInSections[6,2] <-nrow(FurtherCount[FurtherCount[,"Feature"] == "2000Extended",])
CountsInSections[7,2] <-nrow(FurtherCount[FurtherCount[,"Feature"] == "Intergenic",])
CountsInSections[,3] <- (as.numeric(CountsInSections[,2])/sum(as.numeric(CountsInSections[,2]))) * 100

CountsAlongSections <- cbind(midpoints,Counts,Counts/sum(as.numeric(CountsInSections[,2])))
write.table(CountsInSections,file=OutFileCounts,row.names=F,col.names=F,quote=F,sep="\t")
write.table(CountsAlongSections,file=OutFileProfile,row.names=F,col.names=F,quote=F,sep="\t")


