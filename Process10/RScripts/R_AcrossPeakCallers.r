
Bed2GRanges <- function(BedFile,header=F,strand=T){
    if (header){
      StartPos <- grep("Start|start",colnames(BedFile))
      EndPos <- grep("End|end",colnames(BedFile))
      ChrPos <- grep("Chr|chr",colnames(BedFile))
      TempRanges_Bed <- GRanges(seqnames=as.vector(BedFile[-1,ChrPos]),IRanges(start=as.numeric(as.vector(BedFile[-1,StartPos])),end=as.numeric(as.vector(BedFile[-1,EndPos]))),strand=rep("*",nrow(BedFile)-1))

    }else{
      StartPos <- 2
      EndPos <- 3
      ChrPos <- 1
      TempRanges_Bed <- GRanges(seqnames=as.vector(BedFile[,ChrPos]),IRanges(start=as.numeric(as.vector(BedFile[,StartPos])),end=as.numeric(as.vector(BedFile[,EndPos]))),strand=rep("*",nrow(BedFile)))
      TempRanges_Bed

    }
}

Args <- commandArgs(trailingOnly = TRUE)
library(GenomicRanges)

ss <- read.delim("SampleSheet.csv",sep=",")
SampleName <- "JC209"

ss <- read.delim(file.path(getwd(),"SampleSheet.csv"),sep=",")
SampleName <- Args[1]


PeakTables <- vector("list",length=3)
PeakGranges <- vector("list",length=3)

if(length(grep("MacsPeaks",colnames(ss))) > 0){
  if((ss[ss[,"SampleName"] %in% SampleName,"MacsPeaks"] > 0) & !is.na((ss[ss[,"SampleName"] %in% SampleName,"MacsPeaks"]))){
    PeakTables[[1]] <- read.delim(as.vector(ss[ss[,"SampleName"] %in% SampleName,"Macs_name"]),sep="\t",comment.char="#",header=T)
      PeakGranges[[1]] <- Bed2GRanges(PeakTables[[1]],header=T,strand=F)
  }
}
if(length(grep("SicerPeaks",colnames(ss))) > 0){
  if((ss[ss[,"SampleName"] %in% SampleName,"SicerPeaks"] > 0) & !is.na((ss[ss[,"SampleName"] %in% SampleName,"SicerPeaks"]))){
    PeakTables[[2]] <- read.delim(as.vector(ss[ss[,"SampleName"] %in% SampleName,"Sicer_Name"]),sep="\t",header=F)
      PeakGranges[[2]] <- Bed2GRanges(PeakTables[[2]],header=F,strand=F)
  }
}
if(length(grep("TPICsPeaks",colnames(ss))) > 0){
  if((ss[ss[,"SampleName"] %in% SampleName,"TPICsPeaks"] > 0) & !is.na((ss[ss[,"SampleName"] %in% SampleName,"TPICsPeaks"]))){
    PeakTables[[3]] <- read.delim(as.vector(ss[ss[,"SampleName"] %in% SampleName,"TPICS_Name"]),sep=" ",header=F)
      PeakGranges[[3]] <- Bed2GRanges(PeakTables[[3]],header=F,strand=F)
  }
}


TotalMacsPeaks <- length(PeakGranges[[1]])
TotalMacsLength <- sum(width(PeakGranges[[1]]))

TotalTPICsPeaks <- length(PeakGranges[[3]])
TotalTPICsLength <- sum(width(PeakGranges[[3]]))

TotalSicerPeaks  <- length(PeakGranges[[2]])
TotalSicerLength <- sum(width(PeakGranges[[2]]))



MacsVsTPICs <- vector("numeric",length=2)
TotalMacsInTpics <- length(PeakGranges[[1]][PeakGranges[[1]] %in% PeakGranges[[3]]])
TotalTpicsInMacs <- length(PeakGranges[[3]][PeakGranges[[3]] %in% PeakGranges[[1]]])
PercentMacsInTPICs <- (TotalMacsInTpics/TotalMacsPeaks)*100
PercentTPICsInMacs <- (TotalTpicsInMacs/TotalTPICsPeaks)*100
UnionMacsTPICs<- sum(width(reduce(union(PeakGranges[[1]],PeakGranges[[3]]))))
IntersectionMacsTPICS <- sum(width(reduce(intersect(PeakGranges[[1]],PeakGranges[[3]]))))
JacardIndeMacsTPICS <-   IntersectionMacsTPICS/UnionMacsTPICs


MacsVsSicer <- vector("numeric",length=2)
TotalMacsInSicer <- length(PeakGranges[[1]][PeakGranges[[1]] %in% PeakGranges[[2]]])
TotalSicerInMacs <- length(PeakGranges[[2]][PeakGranges[[2]] %in% PeakGranges[[1]]])
PercentMacsInSicer <- (TotalMacsInSicer/TotalMacsPeaks)*100
PercentSicerInMacs <- (TotalSicerInMacs/TotalSicerPeaks)*100
UnionMacsSicer<- sum(width(reduce(union(PeakGranges[[1]],PeakGranges[[2]]))))
IntersectionMacsSicer <- sum(width(reduce(intersect(PeakGranges[[1]],PeakGranges[[2]]))))
JacardIndexMacsSicer <-   IntersectionMacsSicer/UnionMacsSicer


SicerVsTPICs <- vector("numeric",length=2)
TotalTPICsInSicer <- length(PeakGranges[[3]][PeakGranges[[3]] %in% PeakGranges[[2]]])
TotalSicerInTPICs <- length(PeakGranges[[2]][PeakGranges[[2]] %in% PeakGranges[[3]]])
PercentTPICsInSicer <- (TotalTPICsInSicer/TotalTPICsPeaks)*100
PercentSicerInTPICs <- (TotalSicerInTPICs/TotalSicerPeaks)*100
UnionTPICSSicer<- sum(width(reduce(union(PeakGranges[[3]],PeakGranges[[2]]))))
IntersectionTPICSSicer <- sum(width(reduce(intersect(PeakGranges[[3]],PeakGranges[[2]]))))
JacardIndexTPICSSicer <-   IntersectionTPICSSicer/UnionTPICSSicer


Percents <- c(PercentMacsInTPICs,PercentTPICsInMacs,PercentMacsInSicer,PercentSicerInMacs,PercentTPICsInSicer,PercentSicerInTPICs,JacardIndeMacsTPICS,JacardIndexMacsSicer,JacardIndexTPICSSicer)
Names <- c("Macs_In_TPICs","TPICs_In_Macs","Macs_In_Sicer","Sicer_In_Macs","TPICs_In_Sicer","Sicer_In_TPICs","JacardIndex_TPICS_And_Macs","JacardIndex_Sicer_And_Macs","JacardIndex_TPICs_And_Sicer")
TableOut <- cbind(Names,Percents)

write.table(TableOut,file.path(getwd(),"Peaks","PeakProfiles",paste(SampleName,"_AcrossPeakCallers.txt",sep="")),sep="\t")





