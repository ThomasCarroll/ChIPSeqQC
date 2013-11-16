library(ggplot2)
Arguments <- commandArgs(trailingOnly = T)
WkgDir <- Arguments[1]
ss <- read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=",",stringsAsFactors=F)

DispersionProfiles <- dir(file.path(WkgDir,"Coverage"),pattern="*_CoverageDispersion.RData",full.names=T)
ssds <- vector("numeric",length=length(DispersionProfiles))
ginis <- vector("numeric",length=length(DispersionProfiles))
giniAdjusteds <- vector("numeric",length=length(DispersionProfiles))

namesForMat <- vector("character",length=length(DispersionProfiles))
#ReadAmountsForMat <- vector("numeric",length=length(TSSProfilesFiles))
FileNames <- gsub("/.*/","",gsub(".bwa.*","",DispersionProfiles))



for(i in 1:length(DispersionProfiles)){
namesForMat[i]  <- ss[ss[,1] %in% FileNames[i],"SampleName"]
load(DispersionProfiles[i])
ssds[i] <- ssdOfSample
ginis[i] <- giniOfSample[1]
giniAdjusteds[i] <-  giniOfSample[2]
rm(ssdOfSample)
rm(giniOfSample)
}


BigFrame <- cbind(namesForMat,ssds,ginis,giniAdjusteds)
write.table(BigFrame,file.path(WkgDir,"Coverage","SummaryOfCovDispersion.xls"),sep=",")


plot(BigFrame[,2],main="Plot of SSDs of coverage for Sample")
plot(BigFrame[,3],main="Plot of Gini SSDs of coverage for Sample")
plot(BigFrame[,4],main="Plot of adjusted Gini SSDs of coverage for Sample")

#Position <-  1000-5000:0
#BigFrame <- cbind(Position,t(TSSAverageMat))
#temp <- melt(as.data.frame(BigFrame),id.vars=c("Position"))

#p <- ggplot(temp,aes(x=Position,y=value,col=variable))+stat_smooth(se = FALSE)
#save(p,file=file.path(WkgDir,"Coverage","AverageTSSPlot.RData"))
#ggsave(p,file=file.path(WkgDir,"Coverage","AverageTSSPlot.png"))
