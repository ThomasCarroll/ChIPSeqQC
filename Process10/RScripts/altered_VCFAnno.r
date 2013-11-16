PeaksToGenesCRIStyle <- function(ExtGenes,Rit){
require(GenomicRanges)
    MatchedGenes <- ExtGenes[findOverlaps(ExtGenes,Rit)@matchMatrix[,1]]
    MatchedPeaks <- Rit[findOverlaps(ExtGenes,Rit)@matchMatrix[,2]]
	if(length(findOverlaps(ExtGenes,Rit)@matchMatrix[,2]) != 0){ 
	    NotMatchePeaks <- Rit[-findOverlaps(ExtGenes,Rit)@matchMatrix[,2]]
	}else{
		NotMatchePeaks <- Rit
	}
    TempData <- as.data.frame(MatchedPeaks)[-4] ## Removes width part of object which is automatically generated
    colnames(TempData) <- paste("Variant_",colnames(TempData),sep="")
    elementMetadata(MatchedGenes) <- cbind(as.data.frame(elementMetadata(MatchedGenes)),TempData)
    MetaList <- MatchedGenes
    cat("Overlapping to find nearest feature for non-overlapping Peaks")
    #seqnames()
    TempNearestRanges <- GRanges()
    TempNonMatchedByChr <- GRanges()
    for(i in 1:length(unique(seqnames(NotMatchePeaks)))){
      if(any(seqnames(ExtGenes) %in% unique(seqnames(NotMatchePeaks))[i])){
        Index <- nearest(ranges(NotMatchePeaks[seqnames(NotMatchePeaks) %in% unique(seqnames(NotMatchePeaks))[i]]),ranges(ExtGenes[seqnames(ExtGenes) %in% unique(seqnames(NotMatchePeaks))[i]]))
        TempNearGenes <- ExtGenes[seqnames(ExtGenes) %in% unique(seqnames(NotMatchePeaks))[i]][Index]
        TempPeaks <- NotMatchePeaks[seqnames(NotMatchePeaks) %in% unique(seqnames(NotMatchePeaks))[i]]
        TempData2 <- as.data.frame(TempPeaks)[-4] ## Removes width part of object which is automatically generated
        colnames(TempData2) <- paste("Variant_",colnames(TempData2),sep="")
        elementMetadata(TempNearGenes) <- cbind(as.data.frame(elementMetadata(TempNearGenes)),TempData2)
        TempNearestRanges <- c(TempNearestRanges,TempNearGenes)
      }else{
	if(length(NotMatchePeaks) > 0){
        Temp <- NotMatchePeaks[seqnames(NotMatchePeaks) %in% unique(seqnames(NotMatchePeaks))[i]]
        TempNonMatchedByChr <- c(TempNonMatchedByChr,Temp)
	}
      }
      }
    elementMetadata(TempNearestRanges)$Feature <- ("Off_Amplicon")
    levels(elementMetadata(MetaList)$Feature) <-  c(levels(elementMetadata(MetaList)$Feature),"Off_Amplicon")
    MetaList <- c(MetaList,TempNearestRanges)
    Distances <- DistanceTo(MetaList)
    elementMetadata(MetaList) <- cbind(as.data.frame(elementMetadata(MetaList)),Distances)
#    return(MetaList)
    ForPrinting <- as.data.frame(MetaList)
    VarBegging <-  grep("Variant",colnames(ForPrinting))[1]
    Rearranged <- cbind(ForPrinting[,(VarBegging):ncol(ForPrinting)],ForPrinting[,1:(VarBegging-1)])
#    colnames(Rearranged)[grep("Unique_ID",colnames(Rearranged))] <- "Unique_ID"
    ToMergeIn <- as.data.frame(TempNonMatchedByChr)[,-4]
    colnames(ToMergeIn) <- paste("Variant_",colnames(ToMergeIn),sep="")
    TotalVariants <- merge(Rearranged,ToMergeIn,all=T,sort=F)

}

BuildPeakRanges <- function(CRIStyle){
  GRanges(seqnames=seqnames(CRIStyle),ranges=IRanges(start=as.numeric(as.vector(elementMetadata(CRIStyle)$Variant_start)),end=as.numeric(as.vector(elementMetadata(CRIStyle)$Variant_end))))
}



DistanceTo <- function(StrandedFeatures){
  Centredpeaks <- Centred(BuildPeakRanges(StrandedFeatures))
    DistanceTo3PrimeofFeature = end(StrandedFeatures)-start(Centredpeaks)
    DistanceTo5PrimeofFeature = start(StrandedFeatures)-start(Centredpeaks)
    DistanceToCentreOfFeature = start(Centred(StrandedFeatures))-start(Centredpeaks)
  DistancesToFeatures <- cbind(DistanceTo3PrimeofFeature,DistanceTo5PrimeofFeature,DistanceToCentreOfFeature)
  colnames(DistancesToFeatures) <- c("Distance to 3'-end of Feature","Distance to 5'-end of Feature","Distance to Centre of Feature")
  DistancesToFeatures
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


args <- c("GATK.snps.vcf","SLX-4969_amplicons.bed")
args <- c("GATK.indels.vcf","../48_amplicons.bed")


tmp <- readLines(args[1])

startLine <- max(grep("#", tmp))

SNPS <-  read.delim(args[1], skip = startLine-1)

asBed <- gsub(".vcf", ".bed", args[1])
asAnnovar <- gsub(".vcf", ".av", args[1])

##Make the bed file
system(paste("cat", args[1], " | /home/dunnin01/Solexa_Pipeline/pipelines/SV/src/main/vcf2bed.py > ", asBed))

##make annovar input

system(paste("/lustre/reference_data/mib-cri/SV_Pipeline/annovar/convert2annovar.pl -format vcf4", args[1], ">", asAnnovar))

tmpAv <- read.delim(asAnnovar,header=FALSE)
tmpAv <- data.frame(tmpAv, paste("Variant_ID", seq(1:nrow(tmpAv)),sep="_"))
write.table(tmpAv, file=asAnnovar, sep="\t", row.names=FALSE,col.names=FALSE,quote=FALSE)

system(paste("/lustre/reference_data/mib-cri/SV_Pipeline/annovar/annotate_variation.pl -buildver hg19 ",asAnnovar, " /lustre/reference_data/mib-cri/SV_Pipeline/annovar/humandb"))

exonicOut <- paste(asAnnovar, ".exonic_variant_function",sep="")
variantOut <-  paste(asAnnovar, ".variant_function",sep="")

exInf <- read.delim(exonicOut,header=FALSE)
colnames(exInf)[ncol(exInf)] <- "Variant_Unique_ID"

varInf <- read.delim(variantOut,header=FALSE)
colnames(varInf)[ncol(varInf)] <- "Variant_Unique_ID"

SNPpos <- read.delim(asBed, sep="\t",header=FALSE)

SNPmeta <- SNPS[,1:8]
SNPcalls <- SNPS[,10:ncol(SNPS)]


Amps <- read.delim(args[2],sep="\t",header=F)

AmpRanges <- GRanges(seqnames=as.vector(Amps[!Amps[,1]=="",1]),IRanges(Amps[!Amps[,1]=="",2],Amps[!Amps[,1]=="",3]),strand="*",Name=Amps[!Amps[,1]=="",5], Feature = "ON_Amplicon")


SNPRanges<- GRanges(seqnames=as.vector(SNPpos[!SNPpos[,1]=="",1]),IRanges(SNPpos[!SNPpos[,1]=="",2],SNPpos[!SNPpos[,1]=="",3]),strand="*")
if(ncol(SNPS)>3){
  metaData <-  cbind(SNPmeta,paste("Variant_ID",seq(1,nrow(SNPS[!SNPS[,1]=="",])),sep="_"))
  colnames(metaData)[ncol(metaData)] <- "Unique_ID"
  elementMetadata(SNPRanges) <- metaData
}


AnnotatedVars <- PeaksToGenesCRIStyle(AmpRanges,SNPRanges)

annovarOut <- merge(varInf,exInf,all=TRUE,by.x=grep("Variant_Unique_ID", colnames(varInf)), by.y=grep("Variant_Unique_ID", colnames(exInf)),sort=F)

AnnotatedVars <- merge(AnnotatedVars, annovarOut,all=T,sort=F)

write.table(AnnotatedVars, file=gsub(".vcf", ".annotated.txt",args[1]), quote=FALSE,sep="\t",row.names=FALSE)
