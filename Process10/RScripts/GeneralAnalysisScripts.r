### General Analysi Scripts
GetGRanges <- function(LoadFile,simple=F,sepr="\t"){

#Assumes narrow or broad format

    if(class(LoadFile) == "character"){
      RangesTable <- read.delim(LoadFile,sep=sepr,header=T,skip="#")
    }else if(class(LoadFile) == "matrix"){
      RangesTable <- as.data.frame(LoadFile)
    }else{
      RangesTable <- LoadFile
    }

    Chromosomes <- as.vector(RangesTable[,1])
    Start <- as.numeric(as.vector(RangesTable[,2]))
    End <- as.numeric(as.vector(RangesTable[,3]))
    RegionRanges <- GRanges(seqnames=Chromosomes,ranges=IRanges(start=Start,end=End))
    if(simple == F){
    if(ncol(RangesTable) > 4){
      ID <- as.vector(RangesTable[,4])
      Score <- as.vector(RangesTable[,5])
      if(ncol(RangesTable) > 6){
        Strand <- rep("*",nrow(RangesTable))
        RemainderColumn <- as.data.frame(RangesTable[,-c(1:6)])
        elementMetadata(RegionRanges) <- cbind(ID,Score,Strand,RemainderColumn)
      }else{
        elementMetadata(RegionRanges) <- cbind(ID,Score)
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
        print(i)
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


GetSampleAndInputCountsPerPeaks <- function(PeakTable,SS,config="Config"){

  require(Rsamtools)
  require(GenomicRanges)

  GenomicIDsToFind <- gsub("Peak_","",colnames(PeakTable))
  AllSamplesBAMsToCount <- as.vector(SS[SS[,2] %in% GenomicIDsToFind | gsub("-",".",SS[,2]) %in% GenomicIDsToFind,"Processed_bamFileName"])
  AllInputBAMsToCount <- as.vector(unique(SS[SS[,2] %in% as.vector(SS[SS[,2] %in% GenomicIDsToFind | gsub("-",".",SS[,2]) %in% GenomicIDsToFind,"InputToUse"]),"Processed_bamFileName"]))
  AllBamLocations <- file.path(getwd(),"bamFiles",c(AllSamplesBAMsToCount,AllInputBAMsToCount))
  BamFileRes <- BamFileList(AllBamLocations)
  Which <- GRanges(as.vector(PeakTable[,"Peak_seqnames"]),IRanges(start=as.vector(PeakTable[,"Peak_start"]),end=as.vector(PeakTable[,"Peak_end"])))
  AllCounts <- assays(summarizeOverlaps(Which,BamFileRes))$counts
  PeakTableFull <- cbind(PeakTable,AllCounts)
  return(PeakTableFull)
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
  Centre <- GRanges(seqnames=seqnames(RegionRanges),IRanges(start=round((((as.vector(start(ranges(RegionRanges))))+(as.vector(end(ranges(RegionRanges)))))/2)-distance),end=round((((as.vector(start(ranges(RegionRanges))))+(as.vector(end(ranges(RegionRanges)))))/2)+distance)),strand="*",score=elementMetadata(RegionRanges)$Score)
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
    Start <- as.numeric(as.vector(Tommy[,2]))
    End <- as.numeric(as.vector(Tommy[,3]))
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

FullRes <- function(LoadFile,PotentialChromosomes,bwfile,bamfile,toPlot="Flank",distance=100,distanceIn=50,minD=200,OneSidedFlank=T,NameOfFile=NULL,Total=Total){

#Total <- Counts4[5]
#bwfile  <- BigWigList[5]
#bamfile <- BamList[5]
#toPlot="Centre"
#distance=300
#distanceIn=1
#minD=1
#OneSidedFlank=T
#NameOfFile="TempMe"

    require(GenomicRanges)
    require(rtracklayer)
    RegionRangesList <- GetTargetRegion(LoadFile,distance=distance,distanceIn=distanceIn,min=minD,toPlot,PotentialChromosomes)

    RegionRangesList[[1]] <- WithinBWLimitRanges(bamfile,RegionRangesList[[1]])
    RegionRangesList[[1]] <- RemoveOverLapping(RegionRangesList[[1]])
    LoadFile2 <- RegionRangesList[[1]][order(RegionRangesList[[1]])]

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
    if(OneSidedFlank){
     EdgeMeans <- colMeans(rbind(rev(Rightmeans),Leftmeans))
    }else{
      EdgeMeans <- c(Leftmeans,Rightmeans)
    }
    return(EdgeMeans)
    }
    if(toPlot == "Centre"){
          BWSLeft <- BigWigSelection(ranges = RegionRangesList[[1]])
          TempLeft <- import(bwfile,which=RegionRangesList[[1]],selection = BWSLeft)
          CovResLeft <- rep(TempLeft$score,width(TempLeft))
          YouLeft <-   matrix(CovResLeft,ncol=(distance+distance+1),byrow=T)
          LoadFile3 <- YouLeft[order(as.numeric(elementMetadata(LoadFile2)$score),decreasing=T),]
          if(!is.null(NameOfFile)){
            YouLeftT <- (LoadFile3/Total)*1000000
            Temp <- cbind(paste("V_",seq(1:nrow(YouLeftT)),sep=""),YouLeftT)
            colnames(Temp) <- paste("V_",seq(1:ncol(Temp)),sep="")
            write.table(Temp,NameOfFile,sep="\t",row.names=F,quote=F)
            system(paste("/lustre/mib-cri/carrol09/Work/MyPipe/matrix2png-1.2.2/matrix2png -data ",NameOfFile," -mincolor 255:250:250  -range 0:3 > ",paste(NameOfFile,"png",sep="."),sep=""),wait=T,intern=F)
            cat(paste("/lustre/mib-cri/carrol09/Work/MyPipe/matrix2png-1.2.2/matrix2png -data ",NameOfFile," -mincolor 255:250:250  > ",paste(NameOfFile,"png",sep="."),sep=""))
          }
          CentreMeans <- colMeans(YouLeft)
    return(CentreMeans)
    }
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

Centred <- function(GRanges,distance=0){
require(GenomicRanges)
  PeakStart <- start(GRanges)
  PeakWidth <- width(GRanges)
  PeakCentre <-  PeakStart+round(PeakWidth/2)
  start(ranges(GRanges)) <- PeakCentre-distance
  end(ranges(GRanges)) <- PeakCentre+distance
  GRanges
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
    for(i in 1:length(unique(as.vector(seqnames(NotMatchePeaks))))){
	if(any(seqnames(ExtGenes) %in% unique(as.vector(seqnames(NotMatchePeaks)))[i])){
	print(i)
      Index <- nearest(ranges(NotMatchePeaks[seqnames(NotMatchePeaks) %in% unique(as.vector(seqnames(NotMatchePeaks)))[i]]),ranges(ExtGenes[seqnames(ExtGenes) %in% unique(as.vector(seqnames(NotMatchePeaks)))[i]]))
      TempNearGenes <- ExtGenes[seqnames(ExtGenes) %in% unique(as.vector(seqnames(NotMatchePeaks)))[i]][Index]
      TempPeaks <- NotMatchePeaks[seqnames(NotMatchePeaks) %in% unique(as.vector(seqnames(NotMatchePeaks)))[i]]
      TempData2 <- as.data.frame(TempPeaks)[-4] ## Removes width part of object which is automatically generated
      colnames(TempData2) <- paste("Peak",colnames(TempData2),sep="_")
      elementMetadata(TempNearGenes) <- cbind(as.data.frame(elementMetadata(TempNearGenes)),TempData2)
      TempNearestRanges <- c(TempNearestRanges,TempNearGenes)
	}else{
	cat(paste("No Genes found on chromosome -- ",unique(as.vector(seqnames(NotMatchePeaks)))[i],"\n",sep=""))
	}
    }
    elementMetadata(TempNearestRanges)$Feature <- ("Intragenic")
    levels(elementMetadata(MetaList)$Feature) <-  c(levels(elementMetadata(MetaList)$Feature),"Intragenic")
    MetaList <- c(MetaList,TempNearestRanges)
#    elementMetadata(MetaList2)$Feature[length(MetaList):length(MetaList2)] <- rep("Intragenic",length(TempNearestRanges))
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

PeaksToPeaksCRIStyle <- function(ExtGenes,Rit){
require(GenomicRanges)
    MatchedGenes <- ExtGenes[as.matrix(findOverlaps(ExtGenes,Rit))[,1]]
    MatchedPeaks <- Rit[as.matrix(findOverlaps(ExtGenes,Rit))[,2]]
    NotMatchePeaks <- Rit[-as.matrix(findOverlaps(ExtGenes,Rit))[,2]]
    TempData <- as.data.frame(MatchedPeaks)[-4] ## Removes width part of object which is automatically generated
    colnames(TempData) <- paste("Peak",colnames(TempData),sep="_")
    elementMetadata(MatchedGenes) <- cbind(as.data.frame(elementMetadata(MatchedGenes)),TempData)
    MetaList <- MatchedGenes
    DegreesOfOverlap <- as.data.frame(pintersect(ranges(MetaList),IRanges(start=(elementMetadata(MetaList)$Peak_start),end=(elementMetadata(MetaList)$Peak_end)),strand=(elementMetadata(MetaList)$Peak_strand),resolve.empty="start"))
    colnames(DegreesOfOverlap) <- c("OverLap_Start","OverLap_End","BasePairOverLap")
    elementMetadata(MetaList) <- cbind(as.data.frame(elementMetadata(MetaList)),DegreesOfOverlap)
#    elementMetadata(MetaList[elementMetadata(MetaList)$BasePairOverLap == 0])[,c("OverLap_Start","OverLap_End")] <- "NA"
    ##Distance to Stranded beginning and end of feature.
#    elementMetadata(MetaList[elementMetadata(MetaList)$BasePairOverLap == 0])
#    Centre of Peak to Start and End of Gene.... based on strand
#
#    PosFeatures <-  (MetaList[strand(MetaList) == "+"])
#    NegFeatures <-  (MetaList[strand(MetaList) == "-"])
#    PosDistances <- DistanceTo(PosFeatures)
#    elementMetadata(PosFeatures) <- cbind(as.data.frame(elementMetadata(PosFeatures)),PosDistances)
#    NegDistances <- DistanceTo(NegFeatures)
#    elementMetadata(NegFeatures) <- cbind(as.data.frame(elementMetadata(NegFeatures)),NegDistances)
 #   new_MetaList <-  c(PosFeatures,NegFeatures)
    return(MetaList)

}



DistanceTo <- function(StrandedFeatures){
  Centredpeaks <- Centred(BuildPeakRanges(StrandedFeatures))
  if(unique(as.vector(strand(StrandedFeatures))) == "+"){
    DistanceTo3PrimeofFeature = end(StrandedFeatures)-start(Centredpeaks)
    DistanceTo5PrimeofFeature = start(StrandedFeatures)-start(Centredpeaks)
    DistanceToCentreOfFeature = start(Centred(StrandedFeatures))-start(Centredpeaks)
  }
  if(unique(as.vector(strand(StrandedFeatures))) == "-"){
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

BuildOverLapRanges <- function(CRIStyle){
  GRanges(seqnames=seqnames(CRIStyle),ranges=IRanges(start=as.numeric(as.vector(elementMetadata(CRIStyle)$OverLap_Start)),end=as.numeric(as.vector(elementMetadata(CRIStyle)$OverLap_End))))
}

BuildPeakRanges <- function(CRIStyle){
  GRanges(seqnames=seqnames(CRIStyle),ranges=IRanges(start=as.numeric(as.vector(elementMetadata(CRIStyle)$Peak_start)),end=as.numeric(as.vector(elementMetadata(CRIStyle)$Peak_end))))
}


GetPeaksFromSS <- function(SS){
TheSamples <- read.delim(SS,sep=",")
AllMacsPeaks <- GRangesList()
namesForStuff <- vector()
k <- 1
for(i in 1:nrow(TheSamples)){
  FileToRange <- file.path(gsub(".xls",".bed",TheSamples[i,"Macs_name"]))
  if(basename(FileToRange) != "NA"){
    AllMacsPeaks[[k]] <- GetGRanges(FileToRange)
    k <- k +1
    namesForStuff[k] <- as.vector(TheSamples[i,"SampleName"])
    print(i)
  }

}
names(AllMacsPeaks) <- as.vector(namesForStuff)[-1]
return(AllMacsPeaks)
}

RepPerGroup <-  function(RepScoreRanges,SampleSheet,SamplesToMerge){
     SampleSheetJust <- SampleSheet[SampleSheet[,2] %in% SamplesToMerge,c("SampleName","Tissue","Factor","Antibody","Condition_1","Condition_2")]
     MyTable <- paste(SampleSheetJust[,2],SampleSheetJust[,3],SampleSheetJust[,5],SampleSheetJust[,6],sep="_")
     for(k in 1:length(unique(MyTable))){
            NamesToSum <- paste(as.vector(SampleSheetJust[MyTable %in% unique(MyTable)[k],1]))
            TempFrame <- data.frame(rowSums(as.data.frame(elementMetadata(RepScoreRanges)[,colnames(elementMetadata(RepScoreRanges)) %in% gsub("-",".",NamesToSum)])))
            colnames(TempFrame) <- unique(MyTable)[k]
            elementMetadata(RepScoreRanges) <- cbind(as.data.frame(elementMetadata(RepScoreRanges)),TempFrame)
            print(k)
     }
  return(RepScoreRanges)
}

#############
#############

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
  FilteredToGenes <-  Genes_Peaks[!elementMetadata(Genes_Peaks)$Feature %in% "Intragenic"]
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

GetImportantLocations <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  TempDirTemp <- ConfigFile[ConfigFile[,2] %in% "tempdirectory",3]
  WkgDirTemp <- ConfigFile[ConfigFile[,2] %in% "workingdirectory",3]
  LocationsDirTemp <- ConfigFile[ConfigFile[,2] %in% "locationsdirectory",3]
  BamDirTemp <- ConfigFile[ConfigFile[,2] %in% "bamdirectory",3]
  FQDirTemp <- ConfigFile[ConfigFile[,2] %in% "fastqdirectory",3]
  WorkFlowDirTemp <- ConfigFile[ConfigFile[,2] %in% "workflowdirectory",3]
  FragLengthDirTemp <- ConfigFile[ConfigFile[,2] %in% "fraglengthdirectory",3]
  MacsDirTemp <- ConfigFile[ConfigFile[,2] %in% "macsdirectory",3]
  SicerDirTemp <- ConfigFile[ConfigFile[,2] %in% "sicerdirectory",3]
  TPICsDirTemp <- ConfigFile[ConfigFile[,2] %in% "tpicsdirectory",3]
  setClass("ChIPDirLocations", representation(TempDir = "character",WkgDir = "character",LocationsDir = "character",BamDir = "character",FQDir = "character",WorkFlowDir="character",FragLengthDir="character",MacsDir="character",SicerDir="character",TPICsDir="character"))
  PLLocations <- new("ChIPDirLocations",TempDir=TempDirTemp,WkgDir=WkgDirTemp,LocationsDir=LocationsDirTemp,BamDir=BamDirTemp,FQDir=FQDirTemp,WorkFlowDir=WorkFlowDirTemp,FragLengthDir=FragLengthDirTemp,MacsDir=MacsDirTemp,SicerDir=SicerDirTemp,TPICsDir=TPICsDirTemp)
  return(PLLocations)
}


ScorePerGroup <-  function(RepScoreRanges,SampleSheet,SamplesToMerge){
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

RPKMPerGroup <-  function(CountedPeakTable,SampleSheet,SamplesToMerge=NULL){

     UpdatedTable <- CountedPeakTable
     if(!is.null(SamplesToMerge)){
      SampleSheetJust <- SampleSheet[SampleSheet[,2] %in% SamplesToMerge,c("GenomicsID","Tissue","Factor","Antibody","Condition_1","Condition_2")]
     }else{
       SampleSheetJust <- SampleSheet[,c("GenomicsID","Tissue","Factor","Antibody","Condition_1","Condition_2")]
     }
     MyTable <- paste(SampleSheetJust[,3],SampleSheetJust[,5],SampleSheetJust[,6],sep="_")
     TrySumPerGroup <- c("Counts","RPKM")
     for(t in 1:length(TrySumPerGroup)){
       CountedPeakTableTemp <- CountedPeakTable[,grep(paste(TrySumPerGroup[t],"_",sep=""),colnames(CountedPeakTable))]
       colnames(CountedPeakTableTemp) <- gsub(paste(TrySumPerGroup[t],"_",sep=""),"",colnames(CountedPeakTableTemp))
       CleanedColNames <-  gsub("\\.bam$","",colnames(CountedPeakTableTemp))
       CleanedColNames <-  gsub("_Processed$","",CleanedColNames)
       CleanedColNames <-  gsub("\\.bwa.*","",CleanedColNames)

       for(k in 1:length(unique(MyTable))){

              NamesToSum <- paste(as.vector(SampleSheetJust[MyTable %in% unique(MyTable)[k],1]))
              if(length(which(gsub("-",".",CleanedColNames) %in% gsub("-",".",NamesToSum))) > 1){
                TempFrame <- data.frame(rowMeans(CountedPeakTableTemp[,gsub("-",".",CleanedColNames) %in% gsub("-",".",NamesToSum)]))
                colnames(TempFrame) <- paste(TrySumPerGroup[t],"_",unique(MyTable)[k],sep="")
                UpdatedTable <- cbind(UpdatedTable,TempFrame)
                print(k)
              }

       }
     }
  return(UpdatedTable)
}



### Need to add "Counts_" to methods
MakeRPKMPerSample <-  function(CountedPeakTable,SampleSheet,TotalColumn="Unique"){
      Temp <- CountedPeakTable[,grep("Counts_",colnames(CountedPeakTable))]
      LengthOfRegions <- width(IRanges(start=as.vector(PeakTable[,"Peak_start"]),end=as.vector(PeakTable[,"Peak_end"])))
      ##Insert 1 read...
      Temp[Temp == 0] <- 1
      RPK <- Temp/LengthOfRegions
      RPKM <- matrix(nrow=nrow(RPK),ncol=ncol(RPK))
      RPKMOverInput <- matrix(nrow=nrow(RPK),ncol=ncol(RPK))
      colnames(RPK) <- gsub("Counts","RPK",colnames(RPK))
      for(c in 1:length(colnames(RPK))){
        SampleName <- gsub("RPK_","",colnames(RPK)[c])

        SampleName <- gsub("\\.bam$","",SampleName)
        SampleName <- gsub("_Processed$","",SampleName)
        SampleName <- gsub("\\.bwa.*","",SampleName)

        RPKM[,c] <- RPK[,c]/as.numeric(SampleSheet[SampleSheet[,"GenomicsID"]  %in% SampleName,TotalColumn])
      }
      RPKM <- RPKM*1000000
      colnames(RPKM) <- gsub("RPK","RPKM",colnames(RPK))
      CountedPeakTable <-  cbind(CountedPeakTable,RPKM)
}



GetSampleAndInputCountsPerPeaks <- function(PeakTable,SS,config="Config"){

  require(Rsamtools)
  require(GenomicRanges)

  GenomicIDsToFind <- gsub("Peak_","",colnames(PeakTable))
  AllSamplesBAMsToCount <- as.vector(SS[SS[,2] %in% GenomicIDsToFind | gsub("-",".",SS[,2]) %in% GenomicIDsToFind,"Processed_bamFileName"])
  AllInputBAMsToCount <- as.vector(unique(SS[SS[,2] %in% as.vector(SS[SS[,2] %in% GenomicIDsToFind | gsub("-",".",SS[,2]) %in% GenomicIDsToFind,"InputToUse"]),"Processed_bamFileName"]))
  AllBamLocations <- file.path(getwd(),"bamFiles",c(AllSamplesBAMsToCount,AllInputBAMsToCount))
  BamFileRes <- BamFileList(AllBamLocations)
  Which <- GRanges(as.vector(PeakTable[,"Peak_seqnames"]),IRanges(start=as.vector(PeakTable[,"Peak_start"]),end=as.vector(PeakTable[,"Peak_end"])))
  AllCounts <- assays(summarizeOverlaps(Which,BamFileRes))$counts
  PeakTableFull <- cbind(PeakTable,AllCounts)
  return(PeakTableFull)
}



HeatMapAndProfiles <- function(Biggertemp,TempSampleSheet,DistanceToUse=300,PicName,sortBy){

Bigger3 <- cbind(as.vector(Biggertemp[,1]),Biggertemp[,2],Biggertemp[,3],paste(Biggertemp[,1],Biggertemp[,2],Biggertemp[,3],sep="_"),Biggertemp[,4])
colnames(Bigger3) <- c("Peak_seqnames","Peak_start","Peak_end","PeakID","Width")


TempSampleSheet <- TempSampleSheet[order(TempSampleSheet[,"Factor"],TempSampleSheet[,"Tissue"],TempSampleSheet[,"Condition_1"]),]
write.table(TempSampleSheet,"FileOrderForPlotting.csv",sep=",",row.names=F)
AllGroups <- paste(TempSampleSheet[,"Factor"],TempSampleSheet[,"Tissue"],TempSampleSheet[,"Condition_1"],sep=".")
UniqAllGroups <- unique(AllGroups)
Bigger3 <- Bigger3[order(Bigger3[,sortBy],decreasing=T),]

LoadFile <- Bigger3[,c("Peak_seqnames","Peak_start","Peak_end","PeakID",sortBy)]

PotentialChromosomes <- c("chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr23","chr24","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX","chrY")


CountSummit <- 0
LoadFileRanges <- GetGRanges(LoadFile)
    cat("Finding Summits..\n")
   for(k in 1:nrow(TempSampleSheet)){
    cat(paste("Processing files for",TempSampleSheet[k,"SampleName"],sep=" "))
    togetSummits <-   gsub("_peaks.xls","_summits.bed",TempSampleSheet[k,"Macs_name"])
    if(file.exists(togetSummits)){
        CountSummit <- CountSummit+1
        Summits <- GetGRanges(togetSummits)
        FiltTemp <- PeaksToPeaksCRIStyle(LoadFileRanges,Summits)
        FiltTemp <- FiltTemp[order(as.numeric(as.vector(elementMetadata(FiltTemp)$Peak_Score)),decreasing=T)]
        FiltTemp <- FiltTemp[match(unique(elementMetadata(FiltTemp)$ID),elementMetadata(FiltTemp)$ID)]
        SummitOfPeaks <- as.data.frame(cbind(as.vector(elementMetadata(FiltTemp)$ID),round((elementMetadata(FiltTemp)$Peak_start+elementMetadata(FiltTemp)$Peak_end)/2),as.vector(elementMetadata(FiltTemp)$Peak_Score)))
        LoadFile <- merge(LoadFile,SummitOfPeaks,by.x=4,by.y=1,all.x=T,all.y=F)
        if(CountSummit == 1){
        LoadFile <- cbind(LoadFile[,c(2,3,4,1,5:(ncol(LoadFile)-2))],rep("*",nrow(LoadFile)),LoadFile[,c((ncol(LoadFile)-1):(ncol(LoadFile)))])
        }else{
        LoadFile <- cbind(LoadFile[,c(2,3,4,1,5:(ncol(LoadFile)-2))],LoadFile[,c((ncol(LoadFile)-1):(ncol(LoadFile)))])
        }
        colnames(LoadFile)[(ncol(LoadFile)-1):ncol(LoadFile)] <- paste(c("SummitPos","SummitScore"),as.vector(TempSampleSheet[k,"GenomicsID"]),sep="_")
        LoadFileRanges <- GetGRanges(LoadFile)
        print(TempSampleSheet[k,"Macs_name"])
  }

}
SummitsAndAllFrame <- as.data.frame(elementMetadata(LoadFileRanges))
SummitPositionsOnly <- SummitsAndAllFrame[,grep("SummitPos",colnames(SummitsAndAllFrame))]
#data.matrix(SummitPositionsOnly)
MeanSummitCentres <- round(rowMeans(matrix(as.numeric(as.vector(as.matrix(SummitPositionsOnly))),ncol=ncol(SummitPositionsOnly)),na.rm=T))

LoadFileRanges <- GetGRanges(LoadFile)
elementMetadata(LoadFileRanges) <- cbind(as.data.frame(elementMetadata(LoadFileRanges)),MeanSummitCentres)

SummitRanges <- cbind(as.vector(seqnames(LoadFileRanges)),as.data.frame(elementMetadata(LoadFileRanges))$MeanSummitCentres,as.data.frame(elementMetadata(LoadFileRanges))$MeanSummitCentres+1,as.data.frame(elementMetadata(LoadFileRanges))$ID,as.data.frame(elementMetadata(LoadFileRanges))$Score,as.data.frame(elementMetadata(LoadFileRanges))$Strand)


#####Fix this properly
SummitRanges <- SummitRanges[!(SummitRanges[,2]) == "NaN",]

SummitRanges <- SummitRanges[order(as.numeric(as.vector(SummitRanges[,5])),decreasing=T),]

BamList <- as.vector(TempSampleSheet[,"Processed_bamFileName"])
BigWigList <- as.vector(TempSampleSheet[,"BigWig_Files"])



MajorProfile4 <- matrix(nrow=length(BamList),ncol=((DistanceToUse)*2)+1)
MajorProfile44 <- matrix(nrow=length(BamList),ncol=((DistanceToUse)*2)+1)
Counts4 <- vector("numeric",length=length(BamList))
for(i in 1:length(BamList)){
  TempCounts <- sum(as.numeric(matrix(unlist(strsplit(system(paste("/lustre/mib-cri/carrol09/Samtools/samtools-0.1.18/samtools idxstats ",BamList[i],sep=""),intern=T),"\t")),ncol=4,byrow=T)[,3]))
  Counts4[i] <- as.numeric(TempCounts)
  MajorProfile4[i,] <- FullRes(SummitRanges,PotentialChromosomes,BigWigList[i],BamList[i],toPlot="Centre",distance=DistanceToUse,distanceIn=1,min=1,NameOfFile=gsub(".bw",paste(PicName,".Peakcounts",sep=""),BigWigList[i]),Total=Counts4[i])
  MajorProfile44[i,] <- (MajorProfile4[i,]/as.numeric(Counts4[i]))*1000000
  print(i)
}

OutName <- file.path(getwd(),paste(PicName,"AllSamplesHeatmap.png",sep=""))
OutAverage <- file.path(getwd(),paste(PicName,"AveraqeProfile.png",sep=""))
ResizedName <- file.path(dirname(OutAverage),paste("Resized",basename(OutAverage),sep=""))
FullMontage <- file.path(getwd(),paste(PicName,"FullMontage.png",sep=""))


system(paste("/lustre/mib-cri/carrol09/Work/MyPipe/ImageMagick-6.7.7-4/utilities/montage ",paste(paste(gsub(".bw",paste(PicName,".Peakcounts",sep=""),BigWigList),".png",sep=""),sep=" ",collapse=" ")," -mode Concatenate -tile ",length(BigWigList),"x ",OutName,sep=""),intern=F,wait=T)
Dimensions <- unlist(strsplit(unlist(strsplit(system(paste("/lustre/mib-cri/carrol09/Work/MyPipe/ImageMagick-6.7.7-4/utilities/identify ",OutName,sep=""),intern=T)," "))[3],"x"))

SomeColors <- sample(colors(),length(UniqAllGroups))

MaxTemp <- vector("numeric",length(UniqAllGroups))
for(l in 1:length(UniqAllGroups)){
  MaxTemp[l] <- max(colMeans(MajorProfile44[AllGroups %in% UniqAllGroups[l],]))
}
Max <- max(MaxTemp) + (max(MaxTemp)/100)*15

png(OutAverage)
plot(colMeans(MajorProfile44[AllGroups %in% UniqAllGroups[1],]),type="l",col=SomeColors[1],ylab="Normalised_Coverage",ylim=c(0,Max))
for(l in 2:length(UniqAllGroups)){
  lines(colMeans(MajorProfile44[AllGroups %in% UniqAllGroups[l],]),type="l",col=SomeColors[l])
}
legend("topright",legend=UniqAllGroups,fill=SomeColors)
dev.off()



width=as.numeric(Dimensions[2])
height=as.numeric(Dimensions[2])
system(paste("/lustre/mib-cri/carrol09/Work/MyPipe/ImageMagick-6.7.7-4/utilities/convert ",OutAverage," -resize ",width,"x",height,"\\! ",ResizedName,sep=""),intern=F,wait=T)

system(paste("/lustre/mib-cri/carrol09/Work/MyPipe/ImageMagick-6.7.7-4/utilities/montage ",OutName," ",ResizedName," -mode Concatenate -tile 2x ",FullMontage,sep=""),intern=F,wait=T)
}

GetCountsPerSample <- function(PeakTable,SampleSheet){
  GenomicIDsToFind <- gsub("Peak_","",colnames(elementMetadata(PeakTable)))
  AllSamplesBAMsToCount <- as.vector(SampleSheet[SampleSheet[,2] %in% GenomicIDsToFind | gsub("-",".",SampleSheet[,2]) %in% GenomicIDsToFind,"Processed_bamFileName"])
  AllInputBAMsToCount <- as.vector(unique(SampleSheet[SampleSheet[,2] %in% as.vector(SampleSheet[SampleSheet[,2] %in% GenomicIDsToFind | gsub("-",".",SampleSheet[,2]) %in% GenomicIDsToFind,"InputToUse"]),"Processed_bamFileName"]))
  AllBamLocations <- c(AllSamplesBAMsToCount,AllInputBAMsToCount)
  BamFileRes <- BamFileList(AllBamLocations)
  AllCounts <- assays(summarizeOverlaps(PeakTable,BamFileRes))$counts
  return(AllCounts)
}


##############################
##############################

  #MaxineARChIP
setwd("/lustre/mib-cri/carrol09/Work/Theresa_HG18/")

SS <- "/lustre/mib-cri/carrol09/Work/Theresa_HG18/SampleSheetEdited.csv"
SampleSheet <- read.delim(SS,sep=",")


Theresa_AR_Peaks <- GetPeaksFromSS(SS)
Theresa_AR_Peaks_WithRepScore <- CountRepPeaks(Theresa_AR_Peaks)
AR_Peaks_WithRepScore <- RepPerGroup(Theresa_AR_Peaks_WithRepScore,SampleSheet,SampleSheet[,"SampleName"])





TempSampleSheet <- SampleSheet[SampleSheet[,"Condition_1"] %in% c("DHT 10nM","MPA 10nM"),]
Atleast2_453_DHT_Not_In_MPA <-  as.data.frame(AR_Peaks_WithRepScore[elementMetadata(AR_Peaks_WithRepScore)$MDA.MB.453_AR_DHT.10nM_NA > 1 & elementMetadata(AR_Peaks_WithRepScore)$MDA.MB.453_AR_MPA.10nM_NA < 1,])
Atleast2_223_DHT_Not_In_MPA <-  as.data.frame(AR_Peaks_WithRepScore[elementMetadata(AR_Peaks_WithRepScore)$MFM223_AR_DHT.10nM_NA > 1 & elementMetadata(AR_Peaks_WithRepScore)$MFM223_AR_MPA.10nM_NA < 1,])
Atleast2_DHT_223_Not_In_453 <-  as.data.frame(AR_Peaks_WithRepScore[elementMetadata(AR_Peaks_WithRepScore)$MFM223_AR_DHT.10nM_NA > 1 & elementMetadata(AR_Peaks_WithRepScore)$MDA.MB.453_AR_DHT.10nM_NA < 1,])
Atleast2_MPA_223_Not_In_453 <-  as.data.frame(AR_Peaks_WithRepScore[elementMetadata(AR_Peaks_WithRepScore)$MFM223_AR_MPA.10nM_NA > 1 & elementMetadata(AR_Peaks_WithRepScore)$MDA.MB.453_AR_MPA.10nM_NA < 1,])
Atleast2_DHT_453_Not_In_223 <-  as.data.frame(AR_Peaks_WithRepScore[elementMetadata(AR_Peaks_WithRepScore)$MDA.MB.453_AR_DHT.10nM_NA > 1 & elementMetadata(AR_Peaks_WithRepScore)$MFM223_AR_DHT.10nM_NA < 1,])
Atleast2_MPA_453_Not_In_223 <-  as.data.frame(AR_Peaks_WithRepScore[elementMetadata(AR_Peaks_WithRepScore)$MDA.MB.453_AR_MPA.10nM_NA > 1 & elementMetadata(AR_Peaks_WithRepScore)$MFM223_AR_MPA.10nM_NA < 1,])




##Bigger3
##TempSampleSheet

HeatMapAndProfiles(Biggertemp=Atleast2_453_DHT_Not_In_MPA,TempSampleSheet,DistanceToUse=750,PicName="Atleast2_453_DHT_Not_In_MPA",sortBy="Width")
HeatMapAndProfiles(Biggertemp=Atleast2_223_DHT_Not_In_MPA,TempSampleSheet,DistanceToUse=750,PicName="Atleast2_223_DHT_Not_In_MPA",sortBy="Width")
HeatMapAndProfiles(Biggertemp=Atleast2_MPA_223_Not_In_453,TempSampleSheet,DistanceToUse=750,PicName="Atleast2_MPA_223_Not_In_453",sortBy="Width")
HeatMapAndProfiles(Biggertemp=Atleast2_MPA_223_Not_In_453,TempSampleSheet,DistanceToUse=750,PicName="Atleast2_MPA_223_Not_In_453",sortBy="Width")

HeatMapAndProfiles(Biggertemp=Atleast2_DHT_453_Not_In_223,TempSampleSheet,DistanceToUse=750,PicName="Atleast2_DHT_453_Not_In_223",sortBy="Width")
HeatMapAndProfiles(Biggertemp=Atleast2_DHT_453_Not_In_223,TempSampleSheet,DistanceToUse=750,PicName="Atleast2_DHT_453_Not_In_223",sortBy="Width")


######
##############



  #MaxineARChIP
setwd("/lustre/mib-cri/carrol09/Work/NewTheresaSamples/")

SS <- "/lustre/mib-cri/carrol09/Work/NewTheresaSamples/SampleSheet_2.csv"
SampleSheet <- read.delim(SS,sep=",")


Theresa_AR_Peaks <- GetPeaksFromSS(SS)
Theresa_AR_Peaks_WithRepScore <- CountRepPeaks(Theresa_AR_Peaks)
AR_Peaks_WithRepScore <- RepPerGroup(Theresa_AR_Peaks_WithRepScore,SampleSheet,SampleSheet[,"SampleName"])

########
########
##Maybe....
#CountsPerSample <- GetCountsPerSample(AR_Peaks_WithRepScore,SampleSheet)
############

#PeakTable <- AR_Peaks_WithRepScore




############

TempSampleSheet <- SampleSheet[SampleSheet[,"Condition_1"] %in% c("DHT 10nM","MPA 10nM"),]
Atleast2_453_DHT_Not_In_MPA <-  as.data.frame(AR_Peaks_WithRepScore[elementMetadata(AR_Peaks_WithRepScore)$MDA.MB.453_AR_DHT.10nM_NA > 1 & elementMetadata(AR_Peaks_WithRepScore)$MDA.MB.453_AR_MPA.10nM_NA < 1,])
Atleast2_223_DHT_Not_In_MPA <-  as.data.frame(AR_Peaks_WithRepScore[elementMetadata(AR_Peaks_WithRepScore)$MFM223_AR_DHT.10nM_NA > 1 & elementMetadata(AR_Peaks_WithRepScore)$MFM223_AR_MPA.10nM_NA < 1,])
Atleast2_DHT_223_Not_In_453 <-  as.data.frame(AR_Peaks_WithRepScore[elementMetadata(AR_Peaks_WithRepScore)$MFM223_AR_DHT.10nM_NA > 1 & elementMetadata(AR_Peaks_WithRepScore)$MDA.MB.453_AR_DHT.10nM_NA < 1,])
Atleast2_MPA_223_Not_In_453 <-  as.data.frame(AR_Peaks_WithRepScore[elementMetadata(AR_Peaks_WithRepScore)$MFM223_AR_MPA.10nM_NA > 1 & elementMetadata(AR_Peaks_WithRepScore)$MDA.MB.453_AR_MPA.10nM_NA < 1,])
Atleast2_DHT_453_Not_In_223 <-  as.data.frame(AR_Peaks_WithRepScore[elementMetadata(AR_Peaks_WithRepScore)$MDA.MB.453_AR_DHT.10nM_NA > 1 & elementMetadata(AR_Peaks_WithRepScore)$MFM223_AR_DHT.10nM_NA < 1,])
Atleast2_MPA_453_Not_In_223 <-  as.data.frame(AR_Peaks_WithRepScore[elementMetadata(AR_Peaks_WithRepScore)$MDA.MB.453_AR_MPA.10nM_NA > 1 & elementMetadata(AR_Peaks_WithRepScore)$MFM223_AR_MPA.10nM_NA < 1,])




##Bigger3

##TempSampleSheet
write.table(TempSampleSheet,"SampleSheetForHeatmap.csv",sep=",",row.names=F)
write.table(Atleast2_MPA_223_Not_In_453,"Atleast2_MPA_223_Not_In_453.txt",sep="\t",row.names=F)
system(paste("bsub -q bioinformatics -J R_8G -R \"rusage[mem=16000]\"  /lustre/mib-cri/carrol09/Work/MyPipe/R/R-2.15.0/bin/R /lustre/mib-cri/carrol09/Work/NewTheresaSamples/MakeHeatmapScript.r","SampleSheetForHeatmap.csv","Atleast2_MPA_223_Not_In_453","width","250","Atleast2_MPA_223_Not_In_453.txt",sep=" "),wait=F,intern=F)

HeatMapAndProfiles(Biggertemp=Atleast2_453_DHT_Not_In_MPA,TempSampleSheet,DistanceToUse=750,PicName="Atleast2_453_DHT_Not_In_MPA",sortBy="Width")
HeatMapAndProfiles(Biggertemp=Atleast2_223_DHT_Not_In_MPA,TempSampleSheet,DistanceToUse=750,PicName="Atleast2_223_DHT_Not_In_MPA",sortBy="Width")
HeatMapAndProfiles(Biggertemp=Atleast2_MPA_223_Not_In_453,TempSampleSheet,DistanceToUse=750,PicName="Atleast2_MPA_223_Not_In_453",sortBy="Width")
HeatMapAndProfiles(Biggertemp=Atleast2_MPA_223_Not_In_453,TempSampleSheet,DistanceToUse=750,PicName="Atleast2_MPA_223_Not_In_453",sortBy="Width")

HeatMapAndProfiles(Biggertemp=Atleast2_DHT_453_Not_In_223,TempSampleSheet,DistanceToUse=750,PicName="Atleast2_DHT_453_Not_In_223",sortBy="Width")
HeatMapAndProfiles(Biggertemp=Atleast2_DHT_453_Not_In_223,TempSampleSheet,DistanceToUse=750,PicName="Atleast2_DHT_453_Not_In_223",sortBy="Width")


###################################################
####################################################


load("C:\\Users\\carrol09\\Work\\ChIPseqPipeline\\Annotation\\Genes_GRCh37.RData")
Genes <- Genes37
library(BSgenome.Hsapiens.UCSC.hg19)

seqlengths(Genes) <- seqlengths(Hsapiens)[match(names(seqlengths(Genes)),names(seqlengths(Hsapiens)))]


elementMetadata(Genes)$Feature <- "Gene"
GenesAndPromoters <- ExtendGenes(Genes,3000)
seqlengths(GenesAndPromoters) <- seqlengths(Genes)
