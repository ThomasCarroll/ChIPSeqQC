Arguments <- commandArgs(trailingOnly = T)
getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
library(GenomicRanges)
library(ggplot2)


WkgDir <- Arguments[1]
WkgDir <- getwd()
#sampleSheet <- read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=",",stringsAsFactors=F)
OffHistFiles <- dir(path=file.path(WkgDir,"Peaks","PeakProfiles"),pattern="*AllOutSidePeaks.hist$",full.names=T)
namesForGraph <- dir(path=file.path(WkgDir,"Peaks","PeakProfiles"),pattern="*AllOutSidePeaks.hist$",full.names=F)


for(i in 1:length(OffHistFiles)){
   DataIn <- read.delim(OffHistFiles[i],sep="\t")[,c(2,3)]
   print(paste("Merging hist from ",OffHistFiles[i],sep=""))
   if(i == 1){
   BigFrame <- DataIn
   }else{
   BigFrame <- merge(BigFrame,DataIn,by=1,all=T)
   }
}

OnHistFiles <- dir(path=file.path(WkgDir,"Peaks","PeakProfiles"),pattern="*AllInPeaks.hist$",full.names=T)
namesForGraph2 <- dir(path=file.path(WkgDir,"Peaks","PeakProfiles"),pattern="*AllInPeaks.hist$",full.names=F)

for(i in 1:length(OnHistFiles)){
   DataIn <- read.delim(OnHistFiles[i],sep="\t")[,c(2,3)]
   print(paste("Merging hist from ",OnHistFiles[i],sep=""))
   if(i == 1){
   BigFrame2 <- DataIn
   }else{
   BigFrame2 <- merge(BigFrame2,DataIn,by=1,all=T)
   }
}



colnames(BigFrame) <-  c("Depth",namesForGraph)
colnames(BigFrame2) <-  c("Depth",namesForGraph2)


temp <- melt(BigFrame,id="Depth")
temp2 <- melt(BigFrame2,id="Depth")
temp <- cbind(temp,"Outside_peaks")
temp2 <- cbind(temp2,"In_peaks")
colnames(temp)[4] <- "Inside_Outside_Peaks"
colnames(temp2)[4] <- "Inside_Outside_Peaks"

longFrame <- rbind(temp,temp2)
longFrame[,3] <- log2(longFrame[,3])
longFrame[,2] <- gsub("_.*","",longFrame[,2])

ggplot(longFrame,aes(x=Depth,y=value,col=variable,linetype=Inside_Outside_Peaks))+geom_smooth()+ylab("Log2 Base-Pairs")
ggsave(file=file.path(WkgDir,"Peaks","PeakProfiles","OnAndOffPeakCoverage.png"))

ss <- read.delim("SampleSheet.csv",sep=",",header=T,stringsAsFactors=F)
Tempss <- ss[,c("GenomicsID","Filtered")]

OnHistFiles <- dir(path=file.path(WkgDir,"Peaks","PeakProfiles"),pattern="*MergedCounts.bed$",full.names=T)
namesForGraph <- gsub("_MergedCounts.bed","",dir(path=file.path(WkgDir,"Peaks","PeakProfiles"),pattern="*MergedCounts.bed$",full.names=F))

TotalRatio  <- vector("numeric",length=length(OffHistFiles))
Total  <- vector("numeric",length=length(OffHistFiles))
TotalOn  <- vector("numeric",length=length(OffHistFiles))
TotalOff <-  vector("numeric",length=length(OffHistFiles))
TotalOnPercent  <- vector("numeric",length=length(OffHistFiles))
TotalOffPercent <- vector("numeric",length=length(OffHistFiles))

for(i in 1:length(OnHistFiles)){
	TempCount <- read.delim(OnHistFiles[i],sep="\t",header=F)	
	Total[i] <- as.numeric(ss[ss[,c("GenomicsID")] %in% namesForGraph[i],c("Filtered")])
	TotalOn[i] <- sum(TempCount[,4])
	TotalOff[i] <- Total[i]-TotalOn[i]
	TotalRatio[i] <- TotalOn[i]/TotalOff[i]
	TotalOnPercent[i] <- TotalOn[i]/Total[i]
	TotalOffPercent[i] <- TotalOff[i]/Total[i]
}



names(TotalRatio) <- namesForGraph
png(file.path(WkgDir,"Peaks","PeakProfiles","Plot__of__In_Peak__to__Outside_Peak__Ratios.png"))
par(mar=c(10, 4, 4, 4) + 0.1)
plot(TotalRatio,main="Plot of Off_Target to On_Target Ratios",pch=20,cex=4,xlab="",ylab="On_Target/Off_Target Reads")
axis(1,at=1:length(TotalRatio),labels=gsub("_off.*","",names(TotalRatio)),las=2)
dev.off()

png(file.path(WkgDir,"Peaks","PeakProfiles","BarPlot__of__In_Peak__to__Outside_Peak__Percentages.png"))
par(mar=c(10, 4, 4, 4) + 0.1)
barplot(rbind(TotalOnPercent*100,TotalOffPercent*100),beside=F,main="Plot of Off_Target to On_Target Percentages",names.arg=gsub("_off.*","",namesForGraph),xlab="",ylab="Percent Reads",col=c("green","grey"),las=2,ylim=c(0,180),legend=c("On_Target","Off_Target"))
dev.off()

png(file.path(WkgDir,"Peaks","PeakProfiles","BarPlot__of__In_Peak__to__Outside_Peak__Reads.png"))
par(mar=c(10, 8, 4, 4) + 0.1)
barplot(rbind(TotalOn/1000000,TotalOff/1000000),beside=T,main="Plot of Off_Target to On_Target Reads",names.arg=gsub("_off.*","",namesForGraph),xlab="",ylab="Total Reads\n(Millions)",col=c("green","grey"),las=2,legend=c("On_Target","Off_Target"),args.legend = list(x = "topleft"))
dev.off()


PosCountFiles <- dir(path=file.path(WkgDir,"Peaks","PeakProfiles"),pattern="*CountsByPos.bed$",full.names=T)
NegCountFiles <- dir(path=file.path(WkgDir,"Peaks","PeakProfiles"),pattern="*CountsByNeg.bed$",full.names=T)
namesForGraph <- gsub("_CountsByPos.bed","",dir(path=file.path(WkgDir,"Peaks","PeakProfiles"),pattern="*CountsByPos.bed$",full.names=F))

ratios <- vector("list",length=length(namesForGraph))
for(i in 1:length(namesForGraph)){
	PosFile <- PosCountFiles[grep(namesForGraph[i],PosCountFiles)]
	NegFile <- NegCountFiles[grep(namesForGraph[i],NegCountFiles)]
	CountForPos <- read.delim(PosFile,sep="\t",h=F)
	CountForNeg <- read.delim(NegFile,sep="\t",h=F)
	Total <- merge(CountForPos,CountForNeg[,c(4,6)],by.x=4,by.y=1,all=T)
	Total <- cbind(Total,log2(Total[,6]+Total[,7]),log2(Total[,6]/Total[,7]))
	ratios[[i]] <- Total[,9]
	png(file.path(WkgDir,"Peaks","PeakProfiles",paste(namesForGraph[i],"Ratio_Of_PerStrand_CountsVsTotalCounts.png",sep="")))
	plot(Total[,8],Total[,9],ylab="Log2 Ratio Of Reads from Strands",xlab="Log2 Total Reads")
	dev.off()
}
png(file.path(WkgDir,"Peaks","PeakProfiles","BoxplotOfRatios.png"))
boxplot(ratios,names=namesForGraph)
dev.off()

GCFiles <- dir(path=file.path(WkgDir,"Peaks","PeakProfiles"),pattern="*GC.txt$",full.names=T)
CountFiles <- dir(path=file.path(WkgDir,"Peaks","PeakProfiles"),pattern="*_Counts.bed$",full.names=T)
namesForGraph <- gsub("_GC.txt","",dir(path=file.path(WkgDir,"Peaks","PeakProfiles"),pattern="*GC.txt$",full.names=F))
GCScores <- vector("list",length=length(namesForGraph))
for(i in 1:length(GCScores)){
	GCData <- read.delim(GCFiles[i],sep="\t",h=T)
	CountFile <- CountFiles[grep(namesForGraph[i],CountFiles)]
	CountData <- read.delim(CountFile,sep="\t",h=F)
	GCScores[[i]] <- GCData[,7]
	TotalGC <- as.numeric(ss[ss[,c("GenomicsID")] %in% namesForGraph,c("Filtered")])
	GCandCounts <- merge(GCData,CountData[,c(4,6)],by.x=4,by.y=1,all=T)
	GCandCounts <- cbind(GCandCounts,(((((GCandCounts[,16])/(GCandCounts[,14]))*1000)/TotalGC)*1000000))
	png(file.path(WkgDir,"Peaks","PeakProfiles",paste(namesForGraph[i],"RPKM_Of_PeaksVsGCOfPeaks.png",sep="")))
	smoothScatter(GCandCounts[,7]*100,GCandCounts[,17],xlab="GC Content",ylab="RPKM")
	dev.off()	
}
png(file.path(WkgDir,"Peaks","PeakProfiles","BoxplotOfRatios.png"))


boxplot(GCScores,file=file.path(WkgDir,"Peaks","PeakProfiles","GC_Content_Boxplot.png"))

