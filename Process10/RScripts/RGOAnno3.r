# /home/mib-cri/local/bin/Rscript /lustre/mib-cri/carrol09/Work/MyPipe/Process10/RScripts/RGOAnno.r HG18 /lustre/mib-cri/carrol09/Work/PipelinePracticeSet/20100826_RossInnesC_JC_TumorChIPs/Peaks/Macs_Peaks/SLX-1201.250.s_4.bwa.homo_sapiens_Processed_peaks.bed Out.txt Enriched.txt

Args <- c("GRCh37","/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/20111109_RossAdams_DN_HNF1bChIP/Peaks/Macs_Peaks/SLX-4498.739.s_3.bwa.homo_sapiens_Processed_peaks.bed","Out.txt","Enriched.txt")
Args <- c(
"GRCh37",
"/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/Trial/20121114_MontoyaAR_DN_Hes6ChIP/Peaks/Macs/MCF7.bwa_Processed_peaks.bed",
"/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/Trial/20121114_MontoyaAR_DN_Hes6ChIP/Peaks/Macs/MCF7.bwa_Processed_peaks_Annotated.xls",
"/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/Trial/20121114_MontoyaAR_DN_Hes6ChIP/Peaks/Macs/MCF7.bwa_Processed_peaks_Annotated.summary",
"/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/Trial/20121114_MontoyaAR_DN_Hes6ChIP/Peaks/Macs/MCF7.bwa_Processed_peaks_GO_Results.txt"
)
#Sicer = /lustre/mib-cri/carrol09/Work/PipelinePracticeSet/20120626_RobinsonJ_JC_ARchipFoxOE_ProperInput/Peaks/Sicer_Peaks/SLX-5361.875.s_2.bwa.homo_sapiens_Processed/SLX-5361.875.s_2.bwa.homo_sapiens_Processed-W200-G600-islands-summary-FDR.01
#MACS = /lustre/mib-cri/carrol09/Work/PipelinePracticeSet/Trial/20121114_MontoyaAR_DN_Hes6ChIP/Peaks/Macs/MCF7.bwa_Processed_peaks.bed
#TPICS = /lustre/mib-cri/carrol09/Work/PipelinePracticeSet/20120626_RobinsonJ_JC_ARchipFoxOE_ProperInput/Peaks/TPICS_Peaks/SLX-5269.861.s_1.bwa.homo_sapiens_Processed/SLX-5269.861.s_1.bwa.homo_sapiens_Processed_TPICS_Peaks.bed

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


library(GenomicRanges)
print(Args[1])
if(tolower(Args[1]) %in% tolower("HG18")){
  load("/lustre/mib-cri/carrol09/Work/MyPipe/Annotation/Genes_HG18.RData")
  Genes <- Genes36
}

if(tolower(Args[1]) %in% tolower("GRCh37")){
  load("/lustre/mib-cri/carrol09/Work/MyPipe/Annotation/Genes_GRCh37.RData")
  Genes <- Genes37
}

elementMetadata(Genes)$Feature <- "Gene"
GenesAndPromoters <- ExtendGenes(Genes,2000)


Peaks <- read.delim(Args[2],sep="\t",comment.char="#",header=F)
#Peaks <- read.delim("/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/Peaks/Macs_Peaks/SLX-4497.739.s_2.bwa.homo_sapiens_Processed_peaks.bed",sep="\t",comment.char="#")
PeakRanges <- Bed2GRanges(Peaks,header=F,strand=F)
Genes_Peaks <- PeaksToGenesCRIStyle(GenesAndPromoters,PeakRanges)
ToPrint <- as.data.frame(Genes_Peaks)
colnames(ToPrint) <- gsub("seqnames","Chr",colnames(ToPrint))

TotalGenesWithPeaks <- unique(as.vector(ToPrint[!ToPrint[,"Feature"] %in% "Intergenic","SymbolGenes"]))
PeaksInGenes <- ToPrint[!ToPrint[,"Feature"] %in% "Intergenic",]
Peaks <- unique(paste(PeaksInGenes[,9],PeaksInGenes[,10],PeaksInGenes[,11],sep="_"))

Totals <- c(length(TotalGenesWithPeaks),length(Peaks))
write.table(ToPrint,Args[3],sep="\t",row.names=FALSE)
write.table(Totals,Args[4],sep="\t",row.names=FALSE)


FunctionalTestPeaks <- Genes_Peaks[!elementMetadata(Genes_Peaks)$Feature %in% "Intergenic"]

if(length(FunctionalTestPeaks) > 5){
GOTerms <- as.list(org.Hs.egGO2ALLEGS)
Temp <- unlist(as.list(org.Hs.egSYMBOL))
EntrezToSymbol <- cbind(names(Temp),Temp)

EntrezToSymbol <- EntrezToSymbol[EntrezToSymbol[,2] %in% as.vector(elementMetadata(Genes)$SymbolGenes),]
library(GO.db)


TempGO <- as.list(GOTERM)
Terms <- vector("character",length=length(TempGO))
GOIDs <- vector("character",length=length(TempGO))

for(i in 1:length(TempGO)){
  Terms[i] <- Term(TempGO[[i]])
  GOIDs[i] <- GOID(TempGO[[i]])
}
AnnoForGO <- cbind(GOIDs,Terms)

GOSymbols <- vector("list",length=length(GOTerms))
for(i in 1:length(GOTerms)){
  GOSymbols[[i]] <-  unique(EntrezToSymbol[EntrezToSymbol[,1] %in% GOTerms[[i]],2])
  print(i)
}


names(GOSymbols) <- names(GOTerms)
ResMatrix <- matrix(nrow=length(GOTerms),ncol=6)
for(i in 1:length(GOSymbols)){
  print(i)
  LengthOfTerm <- sum(width(reduce(GenesAndPromoters[elementMetadata(GenesAndPromoters)$SymbolGenes %in% GOSymbols[[i]]])))
  NoOfPeaks <- length(PeakRanges)
  NoInTerm <-  length(FunctionalTestPeaks[elementMetadata(FunctionalTestPeaks)$SymbolGenes %in% GOSymbols[[i]]])
  Total <- 2.7e9
  Prior <- LengthOfTerm/Total
  Res <- binom.test(NoInTerm,NoOfPeaks,Prior)$p.value
  ResMatrix[i,] <- c(names(GOTerms)[i],LengthOfTerm,NoInTerm,Prior,NoInTerm/NoOfPeaks,Res)
}
 ResMatrix <- cbind(ResMatrix,p.adjust(ResMatrix[,6],"BH"))
FullMatrix <- merge(AnnoForGO,ResMatrix,by=1,all.x=F,all.y=T)
colnames(FullMatrix) <- c("GOID","Term_Name","Length_Of_Term","Number_Of_Peaks_In_Term","Expected_Proportion_Of_Peaks","Observed_Proportion_Of_Peaks","P_Value","Adjusted_P_Value")
write.table(FullMatrix,Args[5],sep="\t",row.names=FALSE,quote=F)
#write.table(ResMatrix,Args[4],sep="\t",row.names=FALSE)
}else(write.table("",Args[5],sep="",row.names=FALSE,col.names=F,quote=F))

