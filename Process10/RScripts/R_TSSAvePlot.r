library(ggplot2)
Arguments <- commandArgs(trailingOnly = T)
WkgDir <- Arguments[1]
WkgDir <- getwd()
ss <- read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=",",stringsAsFactors=F)
TSSProfilesFiles <- dir(file.path(WkgDir,"Coverage"),pattern="*Processed.RData",full.names=T)
TSSAverageMat <-  matrix(nrow=length(TSSProfilesFiles),ncol=5001)
namesForMat <- vector("character",length=length(TSSProfilesFiles))
ReadAmountsForMat <- vector("numeric",length=length(TSSProfilesFiles))
FileNames <- gsub(".*TSS_AvCov_","",gsub(".bwa.*","",TSSProfilesFiles))
print(FileNames)
print(TSSProfilesFiles)
for(i in 1:length(TSSProfilesFiles)){
namesForMat[i]  <- ss[ss[,1] %in% FileNames[i],"SampleName"]
ReadAmountsForMat[i]  <- ss[ss[,1] %in% FileNames[i],"Filtered"]
load(TSSProfilesFiles[i])
TSSAverageMat[i,] <- TempColMeans#/as.numeric(ReadAmountsForMat[i])
TSSAverageMat[i,] <- TSSAverageMat[i,]/as.numeric(ReadAmountsForMat[i])
rm(TempColMeans)
}
rownames(TSSAverageMat) <- namesForMat 

Position <-  1000-5000:0
BigFrame <- cbind(Position,t(TSSAverageMat))
temp <- melt(as.data.frame(BigFrame),id.vars=c("Position"))

#p <- ggplot(temp,aes(x=Position,y=value,col=variable))+stat_smooth(se = FALSE)
save(p,file=file.path(WkgDir,"Coverage","AverageTSSPlot.RData"))
#ggsave(p,file=file.path(WkgDir,"Coverage","AverageTSSPlot.png"))


for(i in 1:length(unique(temp[,2]))){
png(file.path(getwd(),"Coverage",paste((unique(temp[,2])[i]),".png",sep="")),width=2000,height=750)
  plot(temp[temp[,2] %in% unique(temp[,2])[i],1],temp[temp[,2] %in% unique(temp[,2])[i],3],type="l",xlab="Position relative to TSS",ylab="Normalised Coverage")
dev.off()
}

library(Hmisc)
Maxes <- summarize(temp[,3],by=temp[,2],max)
OrderedMaxes <- Maxes[order(Maxes[,2],decreasing=T),]
Top2 <- as.vector(OrderedMaxes[1:3,1])
OrderedMaxes <- Maxes[order(Maxes[,2],decreasing=F),]
Bottom2 <- as.vector(OrderedMaxes[1:3,1])
ToShow <-  unique(c(Top2,Bottom2))


colToUse <- sample(colours(),length(unique(ToShow)))
png(file.path(getwd(),"Coverage",paste("AverageTSSPlot",".png",sep="")),width=2000,height=750)
plot(temp[,1],temp[,3],type="n",xlab="Position relative to TSS",ylab="Normalised Coverage")
for(i in 1:length(ToShow)){
lines(temp[temp[,2] %in% ToShow[i],1],temp[temp[,2] %in% ToShow[i],3],col=colToUse[i])
}
legend("topleft",legend=ToShow,fill=colToUse)
dev.off()


