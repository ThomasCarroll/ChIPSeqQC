Args <- commandArgs(trailingOnly = TRUE)
library(org.Hs.eg.db)

BuildOverLapRanges <- function(CRIStyle){
  GRanges(seqnames=seqnames(CRIStyle),ranges=IRanges(start=as.numeric(as.vector(elementMetadata(CRIStyle)$OverLap_Start)),end=as.numeric(as.vector(elementMetadata(CRIStyle)$OverLap_End))))
}

BuildPeakRanges <- function(CRIStyle){
  GRanges(seqnames=seqnames(CRIStyle),ranges=IRanges(start=as.numeric(as.vector(elementMetadata(CRIStyle)$Peak_start)),end=as.numeric(as.vector(elementMetadata(CRIStyle)$Peak_end))))
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

GREATTest <- function(Genes_Peaks,SymbolList,UniverseRanges){
  ResList <- list()

  elementMetadata(Genes_Peaks) <-  cbind(as.data.frame(elementMetadata(Genes_Peaks)),paste(elementMetadata(Genes_Peaks)[,c("Peak_seqnames")],elementMetadata(Genes_Peaks)[,c("Peak_start")],sep="_"))
  colnames(elementMetadata(Genes_Peaks))[ncol(elementMetadata(Genes_Peaks))] <- "Peak_IDs"
  FilteredToGenes <-  Genes_Peaks[!elementMetadata(Genes_Peaks)$Feature %in% "Intergenic"]
  TotalPeaks <- length(unique(elementMetadata(Genes_Peaks)$Peak_IDs))
  TotalPeaksInGenes <- length(unique(elementMetadata(FilteredToGenes)$Peak_IDs))
  InList <- matrix(nrow=length(FilteredToGenes),ncol=length(SymbolList))
  TotalGeneLength <- sum(width(reduce(UniverseRanges)))
  TotalLength <- sum(as.numeric(seqlengths(UniverseRanges)))
  for(i in 1:length(SymbolList)){
         RangesOfInterest <- UniverseRanges[elementMetadata(UniverseRanges)$SymbolGenes %in% SymbolList[[i]]]
         PeakRangesInROI <-   FilteredToGenes[elementMetadata(FilteredToGenes)$SymbolGenes %in% SymbolList[[i]]]
         TotalROILength <- sum(width(reduce(RangesOfInterest)))
         TotalPeaksInROI <- length(unique(elementMetadata(PeakRangesInROI)$Peak_IDs))
         PropInRoiGenes <-  TotalPeaksInROI/TotalPeaksInGenes
         PropOfRoiGenes <- TotalROILength/TotalGeneLength

         PropInRoi <-  TotalPeaksInROI/TotalPeaks
         PropOfRoi <- TotalROILength/TotalLength
         
         SigOfGreaterInGenes <- binom.test(TotalPeaksInROI,TotalPeaksInGenes,PropOfRoiGenes,"greater")$p.value
         SigOfGreater <- binom.test(TotalPeaksInROI,TotalPeaks,PropOfRoi,"greater")$p.value
         RatioOfRatio <- PropInRoi/PropOfRoi
         RatioOfRatioGenes <- PropInRoiGenes/PropOfRoiGenes
         ResList[[i]] <- c(TotalROILength,TotalGeneLength,TotalLength,TotalPeaksInROI,TotalPeaksInGenes,TotalPeaks,RatioOfRatio,SigOfGreater,RatioOfRatioGenes,SigOfGreaterInGenes)
         InList[,i] <-  elementMetadata(FilteredToGenes)$SymbolGenes %in% SymbolList[[i]]
  }
  names(ResList) <- names(SymbolList)
  colnames(InList) <-  names(SymbolList)
  elementMetadata(FilteredToGenes) <- cbind(as.data.frame(elementMetadata(FilteredToGenes)),InList)
  AllRes <- list(ResList,FilteredToGenes)
  names(AllRes) <- c("GREAT_Results","Genes_Annotated")
  return(AllRes)
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


## Function For defining TSS regions
GetTSS <- function(FullGeneBounds,distanceLeft,distanceRight=NULL){

if(is.null(distanceRight)){
distanceRight <- distanceLeft
}

  PosGenes <-  (FullGeneBounds[strand(FullGeneBounds) == "+"])
  NegGenes <-  (FullGeneBounds[strand(FullGeneBounds) == "-"])
  NewPosGenes <- GRanges(seqnames=seqnames(PosGenes),IRanges(start=(as.vector(start(ranges(PosGenes))))-distanceLeft,end=(as.vector(start(ranges(PosGenes))))+distanceRight),strand=strand(PosGenes))
  NewNegGenes <- GRanges(seqnames=seqnames(NegGenes),IRanges(start=(as.vector(end(ranges(NegGenes)))-distanceRight),end=(as.vector(end(ranges(NegGenes))))+distanceLeft),strand=strand(NegGenes))
  names(NewPosGenes) <- names(PosGenes)
  names(NewNegGenes) <- names(NegGenes)
  elementMetadata(NewPosGenes) <- elementMetadata(PosGenes)
  elementMetadata(NewNegGenes) <- elementMetadata(NegGenes)
  AllPromoters <- c(NewPosGenes,NewNegGenes)
  return(AllPromoters)
}





GetGRanges <- function(LoadFile,simple=F){

#Assumes narrow or broad format

    if(class(LoadFile) == "character"){
      RangesTable <- read.delim(LoadFile,sep="\t",header=F,skip="#")
    }
    else if(class(LoadFile) == "matrix"){
      RangesTable <- as.data.frame(LoadFile)
    }else{
      RangesTable <- LoadFile
    }

    Chromosomes <- as.vector(RangesTable[,1])
    Start <- as.vector(RangesTable[,2])
    End <- as.vector(RangesTable[,3])
    RegionRanges <- GRanges(seqnames=Chromosomes,ranges=IRanges(start=Start,end=End))
    if(simple == F){
    if(ncol(RangesTable) > 5){
      ID <- as.vector(RangesTable[,4])
      Score <- as.vector(RangesTable[,5])
      Strand <- as.vector(RangesTable[,6])
      
      if(ncol(RangesTable) > 6){
        RemainderColumn <- as.data.frame(RangesTable[,-c(1:6)])
        elementMetadata(RegionRanges) <- cbind(ID,Score,Strand,RemainderColumn)
      }else{
        elementMetadata(RegionRanges) <- cbind(ID,Score,Strand)
      }
    }
    }
    return(RegionRanges)
}


CountRepPeaks <- function(RegionList){
    ReducedRegionList <- reduce(unlist(RegionList))

    for(i in 1:length(RegionList)){
        CountOverLaps <- countGenomicOverlaps(ReducedRegionList,RegionList[[i]])
        elementMetadata(ReducedRegionList) <- cbind(as.data.frame(elementMetadata(ReducedRegionList)),CountOverLaps)
    }
    
    if(!is.null(names(RegionList))){
      colnames(elementMetadata(ReducedRegionList)) <- names(RegionList)
    }
    Temp <- as.data.frame(elementMetadata(ReducedRegionList))
    Temp[Temp != 0] <- 1
    PeakInRep_Num <- rowSums(Temp)
    elementMetadata(ReducedRegionList) <- cbind(Temp,PeakInRep_Num)
    return(ReducedRegionList)
}


GetTargetRegion <- function(LoadFile,distance=100,distanceIn=50,min=200,PotentialChromosomes){
    if(class(LoadFile) == "character"){
      Tommy <- read.delim(LoadFile,sep="\t",header=F,skip="#")
    }else{
      Tommy <- LoadFile
    }
    Tommy <- Tommy[Tommy[,1] %in% PotentialChromosomes,]
    Chromosomes <- as.vector(Tommy[,1])
    Start <- as.vector(Tommy[,2])
    End <- as.vector(Tommy[,3])
    RegionRanges <- GRanges(seqnames=Chromosomes,ranges=IRanges(start=Start,end=End))
    if(ncol(Tommy) == 6){
      IDs <- as.vector(Tommy[,4])
      Strand = as.vector(Tommy[,6])
      Strand[Strand != "+" | Strand != "-"] <- "*"
      Score = as.vector(Tommy[,5])
      strand(RegionRanges) <- Strand
      elementMetadata(RegionRanges)$ID <- IDs
      elementMetadata(RegionRanges)$Score <- Score
    }
    if(ncol(Tommy) > 6){
      IDs <- as.vector(Tommy[,4])
      Strand = as.vector(Tommy[,6])
      Strand[Strand != "+" | Strand != "-"] <- "*"
      Score = as.vector(Tommy[,5])
      extra <- Tommy[,-c(1:6)]
      strand(RegionRanges) <- Strand
      elementMetadata(RegionRanges)$ID <- IDs
      elementMetadata(RegionRanges)$Score <- Score
      elementMetadata(RegionRanges) <- cbind(as.data.frame(elementMetadata(RegionRanges)),extra)
    }
    RegionRangesList <- GetTargetEdges(RegionRanges,distance=distance,distanceIn=distanceIn,min=min)
}



JaccardScore <- function(GRangesA,GRangesB){
    require(GenomicRanges)
     UnionSet <- reduce(c(reduce(GRangesA),reduce(GRangesB)))
     IntersectSet <- reduce(intersect(reduce(GRangesA),reduce(GRangesB)))
     
     LengthUnionSet <- sum(width(UnionSet))
     if(length(IntersectSet) > 0){
      LengthIntersectSet <- sum(width(IntersectSet))
     }else{
      LengthIntersectSet <- 0
     }
     InterOverUnion <- LengthIntersectSet/LengthUnionSet     
     return(InterOverUnion)
}

JaccardOverGRangeslist <- function(GRangeslist){
  ResultMatrix <- matrix(nrow=length(GRangeslist),ncol=length(GRangeslist)+1)
  colnames(ResultMatrix) <- c("Peak_Set",names(GRangeslist))  
  for(i in 1:length(GRangeslist)){
    for(k in 1:length(GRangeslist)){
    print(i)
    print(k)
      ResultMatrix[i,k+1] <- JaccardScore(GRangeslist[[i]],GRangeslist[[k]])    
    }  
  }
  ResultMatrix[,1] <- names(GRangeslist)
  return(ResultMatrix)
}

SimpleOverlapTest <- function(GRangesA,GRangesB,lengths){
      ReducedRangesA <- reduce(GRangesA)
      ReducedRangesB <- reduce(GRangesB)      
      IntersectSet <- reduce(intersect(ReducedRangesA,ReducedRangesB))
      
      ReducedRangesALength <- sum(width(ReducedRangesA))
      ReducedRangesBLength <- sum(width(ReducedRangesB))
      IntersectSetLength <- sum(width(IntersectSet))
      
      
      
      GenomeLength <- sum(as.numeric(lengths))
      Expected <- (sum(width(ReducedRangesA))/GenomeLength) * (sum(width(ReducedRangesB))/GenomeLength)
      
      PValue <-  -pbinom(IntersectSetLength,ReducedRangesALength,Expected,lower.tail = FALSE,log.p=T)
      return(PValue)
}

SimplePointOverlapTest <- function(GRangesA,GRangesB,lengths){
      ReducedRangesA <- reduce(GRangesA)
      ReducedRangesB <- reduce(GRangesB)      
 
      
      ReducedRangesBPoint <- GRanges(seqnames(ReducedRangesB),IRanges(start=round(rowMeans(cbind(start(ReducedRangesB),end(ReducedRangesB)))),end=round(rowMeans(cbind(start(ReducedRangesB),end(ReducedRangesB))))))
      TotalBPointsInA <- length(ReducedRangesBPoint[ReducedRangesBPoint %in% ReducedRangesA])
      TotalBPoints <-  length(ReducedRangesBPoint)
            
      GenomeLength <- sum(as.numeric(lengths))
      Expected <- sum(width(ReducedRangesA))/GenomeLength
      
      PValue <-  -pbinom(TotalBPointsInA,TotalBPoints,Expected,lower.tail = FALSE,log.p=T)
      return(PValue)
}


SimplePointOverlapTestOverGRangeslist <- function(GRangeslist,lengths){
  ResultMatrix <- matrix(nrow=length(GRangeslist),ncol=length(GRangeslist)+1)
  colnames(ResultMatrix) <- c("Peak_Set",names(GRangeslist))  
  for(i in 1:length(GRangeslist)){
    for(k in 1:length(GRangeslist)){
    print(i)
    print(k)
      ResultMatrix[i,k+1] <- SimplePointOverlapTest(GRangeslist[[i]],GRangeslist[[k]],lengths)    
    }  
  }
  ResultMatrix[,1] <- names(GRangeslist)
  return(ResultMatrix)
}


SimpleOverlapTestOverGRangeslist <- function(GRangeslist,lengths){
  ResultMatrix <- matrix(nrow=length(GRangeslist),ncol=length(GRangeslist)+1)
  colnames(ResultMatrix) <- c("Peak_Set",names(GRangeslist))  
  for(i in 1:length(GRangeslist)){
    for(k in 1:length(GRangeslist)){
    print(i)
    print(k)
      ResultMatrix[i,k+1] <- SimpleOverlapTest(GRangeslist[[i]],GRangeslist[[k]],lengths)    
    }  
  }
  ResultMatrix[,1] <- names(GRangeslist)
  return(ResultMatrix)
}

SimpleOverlapTestOverGRangeslist <- function(GRangeslist,lengths){
  ResultMatrix <- matrix(nrow=length(GRangeslist),ncol=length(GRangeslist)+1)
  colnames(ResultMatrix) <- c("Peak_Set",names(GRangeslist))  
  for(i in 1:length(GRangeslist)){
    for(k in 1:length(GRangeslist)){
    print(i)
    print(k)
      ResultMatrix[i,k+1] <- SimpleOverlapTest(GRangeslist[[i]],GRangeslist[[k]],lengths)    
    }  
  }
  ResultMatrix[,1] <- names(GRangeslist)
  return(ResultMatrix)
}


GetTargetEdges <- function(RegionRanges,distance,distanceIn=0,min){
  RegionRanges <- RegionRanges[width(RegionRanges) > distanceIn & width(RegionRanges) > min]
  Left <- GRanges(seqnames=seqnames(RegionRanges),IRanges(start=(as.vector(start(ranges(RegionRanges))))-distance,end=(as.vector(start(ranges(RegionRanges))))+distanceIn),strand="*")
  Right <- GRanges(seqnames=seqnames(RegionRanges),IRanges(start=(as.vector(end(ranges(RegionRanges))))-distanceIn,end=(as.vector(end(ranges(RegionRanges))))+distance),strand="*")
#  Left <- Left[width(Left) > distanceIn & width(Left) > min]
#  Right <- Right[width(Right) > distanceIn & width(Right) > min]
  RegionRangesList <- GRangesList(Left,Right)
}

GetTargetCentres <- function(RegionRanges,distance,min){
  Centre <- GRanges(seqnames=seqnames(RegionRanges),IRanges(start=round((((as.vector(start(ranges(RegionRanges))))+(as.vector(end(ranges(RegionRanges)))))/2)-distance),end=round((((as.vector(start(ranges(RegionRanges))))+(as.vector(end(ranges(RegionRanges)))))/2)+distance)),strand="*")
  RegionRangesList <- GRangesList(Centre)
}


RemoveOverLapping <- function(RegionRanges){
    TempList <- as.matrix(findOverlaps(RegionRanges,RegionRanges))
    if(any(duplicated(TempList[,1]))){
      NewRegionRanges <- RegionRanges[-TempList[,1][duplicated(TempList[,1])]]
    }else{
      NewRegionRanges <- RegionRanges
    }
    return(NewRegionRanges)

}


GetTargetRegion <- function(LoadFile,distance=100,distanceIn=50,min=200,toPlot="Flank",PotentialChromosomes){
    if(class(LoadFile) == "character"){
      Tommy <- read.delim(LoadFile,sep="\t",header=F,skip="#")
    }else{
      Tommy <- LoadFile
    }
    Tommy <- Tommy[Tommy[,1] %in% PotentialChromosomes,]
    Chromosomes <- as.vector(Tommy[,1])
    Start <- as.vector(Tommy[,2])
    End <- as.vector(Tommy[,3])
    RegionRanges <- GRanges(seqnames=Chromosomes,ranges=IRanges(start=Start,end=End))
    if(ncol(Tommy) == 6){
      IDs <- as.vector(Tommy[,4])
      Strand = as.vector(Tommy[,6])
      Strand[Strand != "+" | Strand != "-"] <- "*"
      Score = as.vector(Tommy[,5])
      strand(RegionRanges) <- Strand
      elementMetadata(RegionRanges)$ID <- IDs
      elementMetadata(RegionRanges)$Score <- Score
    }
    if(ncol(Tommy) > 6){
      IDs <- as.vector(Tommy[,4])
      Strand = as.vector(Tommy[,6])
      Strand[Strand != "+" | Strand != "-"] <- "*"
      Score = as.vector(Tommy[,5])
      extra <- Tommy[,-c(1:6)]
      strand(RegionRanges) <- Strand
      elementMetadata(RegionRanges)$ID <- IDs
      elementMetadata(RegionRanges)$Score <- Score
      elementMetadata(RegionRanges) <- cbind(as.data.frame(elementMetadata(RegionRanges)),extra)
    }
    if(toPlot == "Flank"){
      RegionRangesList <- GetTargetEdges(RegionRanges,distance=distance,distanceIn=distanceIn,min=min)
    }
    if(toPlot == "Centre"){
      RegionRangesList <- GetTargetCentres(RegionRanges,distance=distance,min=min)    
    }
    return(RegionRangesList)
}


  WithinBWLimitRanges <- function(bamfile,RegionToFilter){
    require(Rsamtools)
    t <- scanBamHeader(bamfile)[[1]][["targets"]]
    RangeWithin <- GRanges(seqnames=names(t),ranges=IRanges(start=1,end=t),strand="*")
    return(RegionToFilter[findOverlaps(RegionToFilter,RangeWithin,type="within")@queryHits])
  }


#  LoadFile <- "/lustre/mib-cri/dunnin01/data/paired_end/gDNA/Elke/Controls/TruSeq_exome_targeted_regions.hg19.bed.chr"
#  PotentialChromosomes <- c("chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23","chr24","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")
#  bwfile <- "/lustre/mib-cri/dunnin01/data/paired_end/gDNA/Elke/Controls/Data_CRI/BGel.BGEL.952.s_1.bwa.homo_sapiens.bigwig"
#  bamfile <- "/lustre/mib-cri/dunnin01/data/paired_end/gDNA/Elke/Controls/Data_CRI/BGel.BGEL.952.s_1.bwa.homo_sapiens.bam"

FullRes <- function(LoadFile,PotentialChromosomes,bwfile,bamfile,toPlot="Flank",distance=100,distanceIn=50,minD=200){
    require(GenomicRanges)
    require(rtracklayer)
    RegionRangesList <- GetTargetRegion(LoadFile,distance=distance,distanceIn=distanceIn,min=minD,toPlot,PotentialChromosomes)

    RegionRangesList[[1]] <- WithinBWLimitRanges(bamfile,RegionRangesList[[1]])
    RegionRangesList[[1]] <- RemoveOverLapping(RegionRangesList[[1]])

    if(toPlot == "Flank"){
    RegionRangesList[[2]] <- WithinBWLimitRanges(bamfile,RegionRangesList[[2]])
    RegionRangesList[[2]] <- RemoveOverLapping(RegionRangesList[[2]])

    
    BWSLeft <- BigWigSelection(ranges = RegionRangesList[[1]])
    TempLeft <- import(bwfile,which=RegionRangesList[[1]],selection = BWSLeft)
    CovResLeft <- rep(TempLeft$score,width(TempLeft))
    YouLeft <-   matrix(CovResLeft,ncol=(distance+distanceIn+1),byrow=T)
    Leftmeans <- colMeans(YouLeft)

    BWSRight <- BigWigSelection(ranges = RegionRangesList[[2]])
    TempRight <- import(bwfile,which=RegionRangesList[[2]],selection = BWSRight)
    CovResRight <- rep(TempRight$score,width(TempRight))
    YouRight <-   matrix(CovResRight,ncol=(distance+distanceIn+1),byrow=T)
    Rightmeans <- colMeans(YouRight)
    EdgeMeans <- colMeans(rbind(rev(Rightmeans),Leftmeans))
    
    return(EdgeMeans)
    }
    if(toPlot == "Centre"){
          BWSLeft <- BigWigSelection(ranges = RegionRangesList[[1]])
          TempLeft <- import(bwfile,which=RegionRangesList[[1]],selection = BWSLeft)
          CovResLeft <- rep(TempLeft$score,width(TempLeft))
          YouLeft <-   matrix(CovResLeft,ncol=(distance+distance+1),byrow=T)
          CentreMeans <- colMeans(YouLeft)
    return(CentreMeans)    
    }
}


GetPeaksFromSS <- function(SS){
TheSamples <- read.delim(SS,sep=",")
AllMacsPeaks <- GRangesList()
namesForStuff <- vector()
k <- 1
for(i in 1:nrow(TheSamples)){
  FileToRange <- file.path(dirname(SS),"Peaks","Macs",basename(gsub(".xls",".bed",TheSamples[i,"Macs_name"])))
  if(basename(FileToRange) != "NA"){
    AllMacsPeaks[[k]] <- GetGRanges(FileToRange)
    k <- k +1
    namesForStuff[k] <- as.vector(TheSamples[i,"SampleName"])
  }

}
names(AllMacsPeaks) <- as.vector(namesForStuff)[-1]
return(AllMacsPeaks)
}

RepPerGroup <-  function(RepScoreRanges,SampleSheet,SamplesToMerge){
     SampleSheetJust <- SampleSheet[SampleSheet[,1] %in% SamplesToMerge,c("GenomicsID","Tissue","Factor","Antibody","Condition_1","Condition_2")]
     MyTable <- paste(SampleSheetJust[,2],SampleSheetJust[,5],sep="_")
     for(k in 1:length(unique(MyTable))){
            NamesToSum <- paste(as.vector(SampleSheetJust[MyTable %in% unique(MyTable)[k],1]))
            TempFrame <- data.frame(rowSums(as.data.frame(elementMetadata(RepScoreRanges)[,colnames(elementMetadata(RepScoreRanges)) %in% gsub("-",".",NamesToSum)])))
            colnames(TempFrame) <- unique(MyTable)[k]
            elementMetadata(RepScoreRanges) <- cbind(as.data.frame(elementMetadata(RepScoreRanges)),TempFrame)
            print(k)
     }
  return(RepScoreRanges)
}
