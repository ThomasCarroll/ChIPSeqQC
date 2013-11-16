Arguments <- commandArgs(trailingOnly = T)

CovDir <- Arguments[1]

files <- dir(path=CovDir,pattern="*.hist",full.names=T)
forPic <- gsub(".hist","",dir(path=CovDir,pattern="*.hist"))
for (i in 1:length(files)){
  if(i == 1){
      TempIn <- read.delim(files[i],sep="\t",header=F)
      GenomeCov <- TempIn[TempIn[,1] %in% "genome",]
      ToMerge <- cbind(GenomeCov[,2],log10(GenomeCov[,3]))
  }
  if(i > 1){
      TempIn <- read.delim(files[i],sep="\t",header=F)
      GenomeCov <- TempIn[TempIn[,1] %in% "genome",]
      ToMerge2 <- cbind(GenomeCov[,2],log10(GenomeCov[,3]))
      ToMerge <- merge(ToMerge,ToMerge2,by=1,all=T)
  }
}


ColorsForGraph <- sample(colours()[grep("dark",colours())],ncol(ToMerge))
png(file.path(CovDir,"TempPlot.png"))
for(i in 2:ncol(ToMerge)){
if(i == 2){
  plot(ToMerge[,c(1,i)],col=ColorsForGraph[i-1],type="l")
}
if(i > 2){
  lines(ToMerge[,c(1,i)],col=ColorsForGraph[i-1],type="l")
}
}
legend("topright",col=ColorsForGraph,legend=forPic,lty=1)
dev.off()
