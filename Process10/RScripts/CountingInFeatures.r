# /home/mib-cri/local/bin/Rscript /lustre/mib-cri/carrol09/Work/MyPipe/Process10/RScripts/RGOAnno.r HG18 /lustre/mib-cri/carrol09/Work/PipelinePracticeSet/20100826_RossInnesC_JC_TumorChIPs/Peaks/Macs_Peaks/SLX-1201.250.s_4.bwa.homo_sapiens_Processed_peaks.bed Out.txt Enriched.txt

Args <- c("GRCh37","/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/20111109_RossAdams_DN_HNF1bChIP/Peaks/Macs_Peaks/SLX-4498.739.s_3.bwa.homo_sapiens_Processed_peaks.bed","Out.txt","Enriched.txt")




library(Rsamtools)

GetLengthsOfFeatures <- function(GeneSectionsToCount){
AllFeatureTypes <- levels(elementMetadata(GeneSectionsToCount)$Feature)
FeatureTypeLengths <- vector("numeric",length=length(AllFeatureTypes))
for(i in 1:length(AllFeatureTypes)){
    print(i)
    FeatureTypeLengths[i] <- sum(width(reduce(GeneSectionsToCount[elementMetadata(GeneSectionsToCount)$Feature %in% AllFeatureTypes[i]])))
}

FeatureTypeLengths <- c(NA,FeatureTypeLengths)
return(FeatureTypeLengths)
}



#for(l in 1:length(fileList)){
CountInFeatures <- function(File,TestFeatures,lengths){
#  File <- fileList[20]
levels(elementMetadata(GeneSectionsToCount)$Feature)
  library(GenomicRanges)

  print(File)
  #BamBnd <-  readGappedAlignments("/lustre/mib-cri/stark01/20110401_LewisS_AM_HydrxMeth/bamFiles/SLX-3713.607.s_1.bwa.homo_sapiens.bam",which=ForWhich_2,what=c("rname", "strand", "pos", "qwidth"))
  t <- scanBamHeader(File)
  TempLengths <- lengths
  TempLengths <- TempLengths[TempLengths[,1] %in% names(t[[1]]$targets),]
  TempLengths <- TempLengths[!TempLengths[,1] == "chrM",]
  TempLengths[,2] <- TempLengths[,2]-1000
  BigWhich <- GRanges(TempLengths[,1],IRanges(1,TempLengths[,2]))
  print("reading gapped alignments..")
  BamBnd <-  readGappedAlignments(File,param=ScanBamParam(flag=scanBamFlag(isDuplicate=NA),which=BigWhich,what=c("rname", "strand", "pos", "qwidth")))
  print("done")
  BamReads <- granges(BamBnd)
#  widthForReads <- as.numeric(names(sort(table(width(BamBnd)),decreasing=TRUE)[1]))
#  BamReads <- BamReads[width(BamReads) == widthForReads]
  PosTemp <- BamReads[strand(BamReads) == "+"]
  BamPointsPos <- GRanges(seqnames=seqnames(PosTemp),IRanges(start=start(PosTemp),end=start(PosTemp)),strand="+")
  rm(PosTemp)
  NegTemp <- BamReads[strand(BamReads) == "-"]
  BamPointsNeg <- GRanges(seqnames=seqnames(NegTemp),IRanges(start=start(NegTemp),end=start(NegTemp)),strand="-")
  rm(NegTemp)
  rm(BamReads)
  BamPoints <- c(BamPointsPos,BamPointsNeg)
  rm(BamPointsPos)
  rm(BamPointsNeg)
#  TestRes <- vector("numeric",length=length(TotTest))
#  TestRes[1] <- length(BamPoints)
TotTest <- levels(elementMetadata(GeneSectionsToCount)$Feature)
TestRes <- vector("numeric",length=length(TotTest))
  print("Counting...")
  for(i in 1:length(TotTest)){
  print(TotTest[i])
    TestRes[i] <- sum(countOverlaps(TestFeatures[elementMetadata(TestFeatures)$Feature == TotTest[i]],BamPoints))
  }
  TestRes <- c(length(BamPoints),TestRes)
  names(TestRes) <- c("Total",TotTest)
#  ResMatrix[,l] <- TestRes
#  ResMatrix[nrow(ResMatrix),l] <-  length(BamPoints)-sum(TestRes[-1])
#  ResMatrix[1,l] <-  length(BamPoints)
  rm(list = c("BamBnd","BamReads","BamPoints","BamPointsPos","BamPointsNeg"))
  print(TestRes)
#  print(l)
  return(TestRes)
}


Args <- commandArgs(trailingOnly = TRUE)
#library(org.Hs.eg.db)

BuildOverLapRanges <- function(CRIStyle){
  GRanges(seqnames=seqnames(CRIStyle),ranges=IRanges(start=as.numeric(as.vector(elementMetadata(CRIStyle)$OverLap_Start)),end=as.numeric(as.vector(elementMetadata(CRIStyle)$OverLap_End))))
}

BuildPeakRanges <- function(CRIStyle){
  GRanges(seqnames=seqnames(CRIStyle),ranges=IRanges(start=as.numeric(as.vector(elementMetadata(CRIStyle)$Peak_start)),end=as.numeric(as.vector(elementMetadata(CRIStyle)$Peak_end))))
}


ExtendForReadsPeaksIn <- function(Full_Features){
require(GenomicRanges)
  FullGeneBounds <- Full_Features[elementMetadata(Full_Features)$Feature == "Gene"]
  PosGenes <-  (FullGeneBounds[strand(FullGeneBounds) == "+"])
  NegGenes <-  (FullGeneBounds[strand(FullGeneBounds) == "-"])

#  NewPosGenes <- GRanges(seqnames=seqnames(PosGenes),IRanges(start=(as.vector(start(ranges(PosGenes))))-distance,end=(as.vector(end(ranges(PosGenes))))),strand=strand(PosGenes))
  PosGenes <- PosGenes[as.vector(end(ranges(PosGenes)))-((as.vector(start(ranges(PosGenes))))+501) > 0]
 # NewNegGenes <- GRanges(seqnames=seqnames(NegGenes),IRanges(start=(as.vector(start(ranges(NegGenes)))),end=(as.vector(end(ranges(NegGenes))))+distance),strand=strand(NegGenes))

  TSSPos <- GRanges(seqnames=seqnames(PosGenes),IRanges(start=(as.vector(start(ranges(PosGenes))))-500,end=(as.vector(start(ranges(PosGenes)))))+500,strand=strand(PosGenes))
  TSSPos500to2000 <- GRanges(seqnames=seqnames(PosGenes),IRanges(start=(as.vector(start(ranges(PosGenes))))-2000,end=(as.vector(start(ranges(PosGenes)))))-501,strand=strand(PosGenes))
  TSSPos2000to5000 <- GRanges(seqnames=seqnames(PosGenes),IRanges(start=(as.vector(start(ranges(PosGenes))))-5000,end=(as.vector(start(ranges(PosGenes)))))-2001,strand=strand(PosGenes))
  TSSPos5000to10000 <- GRanges(seqnames=seqnames(PosGenes),IRanges(start=(as.vector(start(ranges(PosGenes))))-10000,end=(as.vector(start(ranges(PosGenes))))-5001),strand=strand(PosGenes))
  PosGeneBody <- GRanges(seqnames=seqnames(PosGenes),IRanges(start=(as.vector(start(ranges(PosGenes))))+501,end=(as.vector(end(ranges(PosGenes))))),strand=strand(PosGenes))

  NegGenes <- NegGenes[((as.vector(end(ranges(NegGenes))))-501)-(as.vector(start(ranges(NegGenes)))) > 0]

  TSSNeg <- GRanges(seqnames=seqnames(NegGenes),IRanges(start=(as.vector(end(ranges(NegGenes))))-500,end=(as.vector(end(ranges(NegGenes))))+500),strand=strand(NegGenes))
  TSSNeg500to2000 <- GRanges(seqnames=seqnames(NegGenes),IRanges(start=(as.vector(end(ranges(NegGenes))))+501,end=(as.vector(end(ranges(NegGenes))))+2000),strand=strand(NegGenes))
  TSSNeg2000to5000 <- GRanges(seqnames=seqnames(NegGenes),IRanges(start=(as.vector(end(ranges(NegGenes))))+2001,end=(as.vector(end(ranges(NegGenes))))+5000),strand=strand(NegGenes))
  TSSNeg5000to10000 <- GRanges(seqnames=seqnames(NegGenes),IRanges(start=(as.vector(end(ranges(NegGenes))))+5001,end=(as.vector(end(ranges(NegGenes))))+10000),strand=strand(NegGenes))
 NegGeneBody <- GRanges(seqnames=seqnames(NegGenes),IRanges(start=(as.vector(start(ranges(NegGenes)))),end=(as.vector(end(ranges(NegGenes))))-501),strand=strand(NegGenes))



  #names(NewPosGenes) <- names(PosGenes)
  #names(NewNegGenes) <- names(NegGenes)
  elementMetadata(TSSPos) <- data.frame(Feature=rep("TSS_500Upstream_500Downstream",length(TSSPos)))
  elementMetadata(TSSPos500to2000) <- data.frame(Feature=rep("Promoter_2000Upstream_500Upstream",length(TSSPos500to2000)))
  elementMetadata(TSSPos2000to5000) <- data.frame(Feature=rep("Promoter_5000Upstream_2000Upstream",length(TSSPos2000to5000)))
  elementMetadata(TSSPos5000to10000) <- data.frame(Feature=rep("Promoter_10000Upstream_5000Upstream",length(TSSPos5000to10000)))
  elementMetadata(PosGeneBody) <-  data.frame(Feature=rep("GeneBody_minus_TSS",length(PosGeneBody)))

  elementMetadata(TSSNeg) <- data.frame(Feature=rep("TSS_500Upstream_500Downstream",length(TSSNeg)))
  elementMetadata(TSSNeg500to2000) <- data.frame(Feature=rep("Promoter_2000Upstream_500Upstream",length(TSSNeg500to2000)))
  elementMetadata(TSSNeg2000to5000) <- data.frame(Feature=rep("Promoter_5000Upstream_2000Upstream",length(TSSNeg2000to5000)))
  elementMetadata(TSSNeg5000to10000) <- data.frame(Feature=rep("Promoter_10000Upstream_5000Upstream",length(TSSNeg5000to10000)))
  elementMetadata(NegGeneBody) <-  data.frame(Feature=rep("GeneBody_minus_TSS",length(NegGeneBody)))


  AllGenesInSections <- c(TSSPos,TSSPos500to2000,TSSPos2000to5000,TSSPos5000to10000,PosGeneBody,TSSNeg,TSSNeg500to2000,TSSNeg2000to5000,TSSNeg5000to10000,NegGeneBody)

 return(AllGenesInSections)
}



ExtendGenes <- function(Full_Features,distance){
require(GenomicRanges)
  FullGeneBounds <- Full_Features[elementMetadata(Full_Features)$Feature == "Gene"]
  PosGenes <-  (FullGeneBounds[strand(FullGeneBounds) == "+"])
  NegGenes <-  (FullGeneBounds[strand(FullGeneBounds) == "-"])
  NewPosGenes <- GRanges(seqnames=seqnames(PosGenes),IRanges(start=(as.vector(start(ranges(PosGenes))))-distance,end=(as.vector(end(ranges(PosGenes))))),strand=strand(PosGenes))
  NewNegGenes <- GRanges(seqnames=seqnames(NegGenes),IRanges(start=(as.vector(start(ranges(NegGenes)))),end=(as.vector(end(ranges(NegGenes))))+distance),strand=strand(NegGenes))
  #names(NewPosGenes) <- names(PosGenes)
  #names(NewNegGenes) <- names(NegGenes)
  elementMetadata(NewPosGenes) <- elementMetadata(PosGenes)
  elementMetadata(NewNegGenes) <- elementMetadata(NegGenes)
  AllPromoters <- c(NewPosGenes,NewNegGenes)
  elementMetadata(AllPromoters)$Feature <- paste(distance,"Extended",sep="")
  #return(AllPromoters)
  New_Full_Features <- c(AllPromoters)
#  New_Full_Features <- New_Full_Features[order(ranges(New_Full_Features))]
  New_Full_Features
}




PeaksToGenesCRIStyle <- function(ExtGenes,Rit){
require(GenomicRanges)
    MatchedGenes <- ExtGenes[as.matrix(findOverlaps(ExtGenes,Rit))[,1]]
    MatchedPeaks <- Rit[as.matrix(findOverlaps(ExtGenes,Rit))[,2]]
    NotMatchePeaks <- Rit[-as.matrix(findOverlaps(ExtGenes,Rit))[,2]]
    TempData <- as.data.frame(MatchedPeaks)[-4] ## Removes width part of object which is automatically generated
    colnames(TempData) <- paste("Peak",colnames(TempData),sep="_")
    elementMetadata(MatchedGenes) <- cbind(as.data.frame(elementMetadata(MatchedGenes)),TempData)
    MetaList <- MatchedGenes
    cat("Overlapping to find nearest feature for non-overlapping Peaks")
    #seqnames()
    TempNearestRanges <- GRanges()
    for(i in 1:length(unique(seqnames(NotMatchePeaks)))){
	if(any(seqnames(ExtGenes) %in% unique(seqnames(NotMatchePeaks))[i])){
	print(i)
      Index <- nearest(ranges(NotMatchePeaks[seqnames(NotMatchePeaks) %in% unique(seqnames(NotMatchePeaks))[i]]),ranges(ExtGenes[seqnames(ExtGenes) %in% unique(seqnames(NotMatchePeaks))[i]]))
      TempNearGenes <- ExtGenes[seqnames(ExtGenes) %in% unique(seqnames(NotMatchePeaks))[i]][Index]
      TempPeaks <- NotMatchePeaks[seqnames(NotMatchePeaks) %in% unique(seqnames(NotMatchePeaks))[i]]
      TempData2 <- as.data.frame(TempPeaks)[-4] ## Removes width part of object which is automatically generated
      colnames(TempData2) <- paste("Peak",colnames(TempData2),sep="_")
      elementMetadata(TempNearGenes) <- cbind(as.data.frame(elementMetadata(TempNearGenes)),TempData2)
      TempNearestRanges <- c(TempNearestRanges,TempNearGenes)
	}else{
	cat(paste("No Genes found on chromosome -- ",unique(seqnames(NotMatchePeaks))[i],"\n",sep=""))
	}
    }
    elementMetadata(TempNearestRanges)$Feature <- ("Intergenic")
    levels(elementMetadata(MetaList)$Feature) <-  c(levels(elementMetadata(MetaList)$Feature),"Intergenic")
    MetaList <- c(MetaList,TempNearestRanges)
#    elementMetadata(MetaList2)$Feature[length(MetaList):length(MetaList2)] <- rep("Intergenic",length(TempNearestRanges))
    DegreesOfOverlap <- as.data.frame(pintersect(ranges(MetaList),IRanges(start=(elementMetadata(MetaList)$Peak_start),end=(elementMetadata(MetaList)$Peak_end)),strand=(elementMetadata(MetaList)$Peak_strand),resolve.empty="start"))
    colnames(DegreesOfOverlap) <- c("OverLap_Start","OverLap_End","BasePairOverLap")
    elementMetadata(MetaList) <- cbind(as.data.frame(elementMetadata(MetaList)),DegreesOfOverlap)
    elementMetadata(MetaList[elementMetadata(MetaList)$BasePairOverLap == 0])[,c("OverLap_Start","OverLap_End")] <- "NA"
    ##Distance to Stranded beginning and end of feature.
#    elementMetadata(MetaList[elementMetadata(MetaList)$BasePairOverLap == 0])
#    Centre of Peak to Start and End of Gene.... based on strand
#
    PosFeatures <-  (MetaList[strand(MetaList) == "+"])
    NegFeatures <-  (MetaList[strand(MetaList) == "-"])
    PosDistances <- DistanceTo(PosFeatures)
    elementMetadata(PosFeatures) <- cbind(as.data.frame(elementMetadata(PosFeatures)),PosDistances)
    NegDistances <- DistanceTo(NegFeatures)
    elementMetadata(NegFeatures) <- cbind(as.data.frame(elementMetadata(NegFeatures)),NegDistances)
    new_MetaList <-  c(PosFeatures,NegFeatures)
    return(new_MetaList)

}

DistanceTo <- function(StrandedFeatures){
  Centredpeaks <- Centred(BuildPeakRanges(StrandedFeatures))
  if(unique(strand(StrandedFeatures)) == "+"){
    DistanceTo3PrimeofFeature = end(StrandedFeatures)-start(Centredpeaks)
    DistanceTo5PrimeofFeature = start(StrandedFeatures)-start(Centredpeaks)
    DistanceToCentreOfFeature = start(Centred(StrandedFeatures))-start(Centredpeaks)
  }
  if(unique(strand(StrandedFeatures)) == "-"){
    DistanceTo3PrimeofFeature = start(StrandedFeatures)-start(Centredpeaks)
    DistanceTo5PrimeofFeature = end(StrandedFeatures)-start(Centredpeaks)
    DistanceToCentreOfFeature = start(Centred(StrandedFeatures))-start(Centredpeaks)
  }
  DistancesToFeatures <- cbind(DistanceTo3PrimeofFeature,DistanceTo5PrimeofFeature,DistanceToCentreOfFeature)
  colnames(DistancesToFeatures) <- c("Distance to 3'-end of Feature","Distance to 5'-end of Feature","Distance to Centre of Feature")
  DistancesToFeatures
}

GetPromoters <- function(FullGeneBounds,distance){
require(GenomicRanges)
  PosGenes <-  (FullGeneBounds[strand(FullGeneBounds) == "+"])
  NegGenes <-  (FullGeneBounds[strand(FullGeneBounds) == "-"])
  NewPosGenes <- GRanges(seqnames=seqnames(PosGenes),IRanges(start=(as.vector(start(ranges(PosGenes))))-distance,end=(as.vector(start(ranges(PosGenes))))-1),strand=strand(PosGenes))
  NewNegGenes <- GRanges(seqnames=seqnames(NegGenes),IRanges(start=(as.vector(end(ranges(NegGenes)))+1),end=(as.vector(end(ranges(NegGenes))))+distance),strand=strand(NegGenes))
  names(NewPosGenes) <- names(PosGenes)
  names(NewNegGenes) <- names(NegGenes)
  AllPromoters <- c(NewPosGenes,NewNegGenes)
  return(AllPromoters)
}



Bed2GRanges <- function(BedFile,header=F,strand=T){
    if (header){
      StartPos <- grep("Start|start",colnames(BedFile))
      EndPos <- grep("End|end",colnames(BedFile))
      ChrPos <- grep("Chr|chr",colnames(BedFile))
    }else{
      StartPos <- 2
      EndPos <- 3
      ChrPos <- 1
    }
      TempRanges_Bed <- GRanges(seqnames=as.vector(BedFile[,ChrPos]),IRanges(start=as.numeric(as.vector(BedFile[,StartPos])),end=as.numeric(as.vector(BedFile[,EndPos]))),strand=rep("*",nrow(BedFile)))
      if(ncol(BedFile) > 3){
		elementMetadata(TempRanges_Bed) <- BedFile[,-c(ChrPos,StartPos,EndPos)]
	}
      TempRanges_Bed
}

Centred <- function(GRanges,distance=0){
require(GenomicRanges)
  PeakStart <- start(GRanges)
  PeakWidth <- width(GRanges)
  PeakCentre <-  PeakStart+round(PeakWidth/2)
  start(ranges(GRanges)) <- PeakCentre-distance
  end(ranges(GRanges)) <- PeakCentre+distance
  GRanges
}


Bed2GeneGRanges <- function(BedFile,header=F,strand=T){
    if (header){
      StartPos <- grep("Start|start",colnames(BedFile))
      EndPos <- grep("End|end",colnames(BedFile))
      ChrPos <- grep("Chr|chr",colnames(BedFile))
      StrandPos <- grep("Strand|strand",colnames(BedFile))
      
    }else{
      StartPos <- 2
      EndPos <- 3
      ChrPos <- 1
      StrandPos <- 4
    } 
      TempRanges_Bed <- GRanges(seqnames=as.vector(BedFile[,ChrPos]),IRanges(start=as.numeric(as.vector(BedFile[,StartPos])),end=as.numeric(as.vector(BedFile[,EndPos]))),strand=as.vector(BedFile[,StrandPos]))
		    elementMetadata(TempRanges_Bed) <- BedFile[,-c(ChrPos,StartPos,EndPos,StrandPos)] 
        colnames(elementMetadata(TempRanges_Bed))[1] <- "SymbolGenes" 
      TempRanges_Bed
}


####################################################################################
####################################################################################

library(GenomicRanges)
print(Args[1])

File <- Args[2]
WkgDir <- Args[3]

  if(any(grepl(".RData",Args[1]))){
    Genes <- get(load(Args[1]))
  }
  if(any(grepl(".bed",Args[1]))){
    GenesTemp <- read.delim(Args[1],sep="\t",comment.char="#",header=T)
    Genes <- Bed2GeneGRanges(GenesTemp,header=T,strand=T)
  }



print("Here1")

#GeneLengths <- read.delim("/lustre/mib-cri/carrol09/MyPipe/bedFiles/hg19.txt",sep="\t",header=F)
GeneLengths <- read.delim(Args[4],sep="\t",header=F)

print("Here")

elementMetadata(Genes)$Feature <- "Gene"
GenesAndPromoters <- ExtendGenes(Genes,10000)
print("Here5")
TempRepeats <- table(elementMetadata(GenesAndPromoters[as.matrix(findOverlaps(GenesAndPromoters,GenesAndPromoters))[,1],])$SymbolGenes)
CrossLappingGenes <- names(TempRepeats[TempRepeats > 2])
NonLappingGenes <-  Genes[!elementMetadata(Genes)$SymbolGenes %in% CrossLappingGenes,]

print("Here2")
GeneSectionsToCount <- ExtendForReadsPeaksIn(NonLappingGenes)

print("Here")

TotTest <- c("Total",levels(elementMetadata(GeneSectionsToCount)$Feature),"Intergenic")

ResMatrix <- matrix(nrow=length(TotTest),ncol=1)

Counts <- CountInFeatures(File,GeneSectionsToCount,GeneLengths)
FeatureTypeLengths <- GetLengthsOfFeatures(GeneSectionsToCount)

ForExport <- rbind(Counts,FeatureTypeLengths)
write.table(ForExport,file=file.path(WkgDir,gsub("/.*/","",gsub("\\.bam","_ReadCountsInFeatures.txt",File))),sep="\t")
