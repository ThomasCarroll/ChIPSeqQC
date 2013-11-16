### Load in Libraries
library(parallel)
library(GenomicRanges)
#library(BSgenome.Hsapiens.UCSC.hg18)
library(Rsamtools)

CoverageFeature <- function(fileList,WindowSizeLeft,WindowSizeRight=NULL,prefix="",selectedGenes,WkgDir){
#require(Rlsf)
#require(multicore)
require(GenomicRanges)
#require(BSgenome.Hsapiens.UCSC.hg18)
require(Rsamtools)


##
if(is.null(WindowSizeRight)){
WindowSizeRight <- WindowSizeLeft
cat("No right window size set..using equal distances around feature.\n")
}

ReadLength <- 100
cat("Set read length to 100bp for windows around features.\n")


cat("Loading in annotation....\n")

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


## Filter Features [Genes/Exons/Introns etc] to just Gene positions
FullGenes <- selectedGenes
## Define start points...


cat("Defining TSSs....\n")

TSS <- GetTSS(FullGenes,WindowSizeLeft+ReadLength,WindowSizeRight+ReadLength)
TSS_2 <- TSS
## Rename Chromosome names [eg. 1 to chr1]
#TSS_2 <- GRanges(seqnames=paste("chr",as.vector(seqnames(TSS)),sep=""),ranges(TSS),strand=strand(TSS))


##Scan in File using
Params <- ScanBamParam(which=TSS_2,what=c("rname", "strand", "pos", "qwidth"))
TSS <- GetTSS(FullGenes,0)
TSS_3 <- TSS
elementMetadata(TSS_3) <- elementMetadata(TSS)
cat("...Done\n")
## Read in a list of Santiago's files
rm(selectedGenes)
#rm(TSS_2)
rm(TSS)
#####
#if(!is.null(selectedGenes)){
#cat("Trimming down gene list to those specified...\n")
#c#at(paste("",length(selectedGenes)," genes provided\n",sep=""))
#  TSS_3 <- TSS_3[elementMetadata(TSS_3)$Gene_Symbol %in% selectedGenes]
#cat(paste("",length(TSS_3)," genes selected\n",sep=""))
#}

GetCoverage <- function(file,TSS_3,direct,prefix){
require(Rlsf)
#require(multicore)
require(GenomicRanges)
#require(BSgenome.Hsapiens.UCSC.hg18)
require(Rsamtools)


  deltaConvolve <- function(Chromosome,TSS_3,TSS_2,ReadRanges, left, right){
require(Rlsf)
#require(multicore)
require(GenomicRanges)
#require(BSgenome.Hsapiens.UCSC.hg18)
require(Rsamtools)
  ######
##  Need to make this strand dependent too ..as in left and right are no longer symmetrical

    print(as.vector(Chromosome))
    TempTSS_3 <- TSS_3[seqnames(TSS_3) %in% as.vector(Chromosome)]
    TempReadRanges <- ReadRanges[seqnames(ReadRanges) %in% as.vector(Chromosome)]
    TempCoverage <- coverage(ranges(TempReadRanges),width=seqlengths(TempReadRanges)[names(seqlengths(TempReadRanges)) %in% as.vector(Chromosome)])
    winCentres <- as.vector(start(TempTSS_3))
    names(winCentres) <- elementMetadata(TempTSS_3)$Gene_Symbol
    strand <- as.vector(strand(TempTSS_3))
    result <- rep(0, right + left + 1)
    count <- 0
    for( i in 1:length(winCentres) ) {
    count <- count+1
    if(strand[i] == "-"){
    leftTemp <- right
    rightTemp <- left
    }
    else{
    leftTemp <- left
    rightTemp <- right
    }


#    print(winCentres[i])
#    print(winCentres[i] + left)
    v <- as.vector(window(TempCoverage,winCentres[i] - leftTemp,
    winCentres[i] + rightTemp))
    if(strand[i] == "-"){v <- rev(v)}
    #v <- rev(v)
#    print(paste("MinorLoop",i,sep=""))
    result <- rbind(result,v)
    }
  cat("Finished on Chromosome..",as.vector(Chromosome),"..processed ",count," features from a supplied",length(winCentres),"and produced",nrow(result)," results\n")
#  write.table(paste("Finished on Chromosome..",as.vector(Chromosome),"..processed ",count," features from a supplied ",length(winCentres)," features\n",sep=""),file=paste(direct,"/TempProcessedFile.txt",sep=""),append=T)
  result <- result[-1,]
  rownames(result) <- names(winCentres)
  result
  }



  cat("Reading in file......\n")
if(file.exists(paste(file,".bai",sep=""))){
  BamBnd <-  readGappedAlignments(file,param=ScanBamParam(which=TSS_2,what=c("rname", "strand", "pos", "qwidth")))
if(length(BamBnd) > 1000){
  cat("...Done\n")
  #BamBnd <-  readGappedAlignments("/lustre/mib-cri/carrol09/Work/ToBreak/20110526_TheodorouV_JC_MYCE2F1/bamFiles/SLX-2574.433.s_2.bwa.homo_sapiens_exclude_f_JC133_f.bam",which=TSS_2,what=c("rname", "strand", "pos", "qwidth"))
  ReadRanges <- granges(BamBnd)
#  Temp <- vector("list",length=length(seqlevels(TSS_3)))
#  for(i in 1:length(seqlevels(TSS_3))){
#     TempTSS_3 <- TSS_3[seqnames(TSS_3) %in% unique(seqnames(TSS_3))[i]]
#     TempReadRanges <- ReadRanges[seqnames(ReadRanges) %in% unique(seqnames(TSS_3))[i]]
#     TempCoverage <- coverage(ranges(TempReadRanges),width=seqlengths(TempReadRanges)[names(seqlengths(TempReadRanges)) %in% unique(seqnames(TSS_3))[i]])
     #Temp[[i]] <- deltaConvolve(TempCoverage,as.vector(start(TempTSS_3)),as.vector(strand(TempTSS_3)),-WindowSize,WindowSize)
     Temp <- mclapply(unique(seqnames(TSS_3))[1:4],deltaConvolve,TSS_3=TSS_3,TSS_2=TSS_2,ReadRanges,WindowSizeLeft,WindowSizeRight,mc.cores=4)
#     print(i)
  #   return(Temp)
#  }

  Temp2 <- Temp[[1]]
  print("1")
  for(i in 2:length(Temp)){
  if(!is.null(ncol(Temp[[i]]))){
  if(ncol(Temp[[i]]) > 1){
    Temp2 <- rbind(Temp2,Temp[[i]])
    }
  }
    print(i)
  }
  print(Temp2)
  TempColMeans <- colMeans(Temp2)
  save(TempColMeans,file=paste(direct,"/",prefix,gsub("/.*/","",gsub(".bam","",file)),".RData",sep=""))
}else{
TempColMeans <- rep(0,(WindowSizeLeft+WindowSizeRight+1))
save(TempColMeans,file=paste(direct,"/",prefix,gsub("/.*/","",gsub(".bam","",file)),".RData",sep=""))
}
}else{
TempColMeans <- rep(0,(WindowSizeLeft+WindowSizeRight+1))
save(TempColMeans,file=paste(direct,"/",prefix,gsub("/.*/","",gsub(".bam","",file)),".RData",sep=""))
}
}
  cat("Running Coverage...\n")
  GetCoverage(fileList,TSS_3=TSS_3,WkgDir,prefix)
  cat("..Done\n")
}


Arguments <- commandArgs(trailingOnly = T)
#Arguments <- c("/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/AnotherVerison/20111109_RossAdams_DN_HNF1bChIP/bamFiles/SLX-4500.739.s_5.bwa.homo_sapiens_Processed.bam","/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/AnotherVerison/20111109_RossAdams_DN_HNF1bChIP/","GRCh37")

getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))


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


OnlyGenesInLims <- function(Genes,ChrLims){
  TrimmedGenes <- vector("list",length=length(ChrLims))
  ChrNames <-  as.character(as.vector(ChrLims[,1]))
  ChrLengths <- as.numeric(as.vector(ChrLims[,2]))
  for(i in 1:length(unique(ChrNames))){
      ChrGenes  <- Genes[seqnames(Genes) %in% ChrNames[i]]
	if(length(ChrGenes[start(ChrGenes)-4000 > 0 & end(ChrGenes)+1000 < ChrLengths[i]]) > 0){
   	   TrimmedGenes[[i]] <- ChrGenes[start(ChrGenes)-4000 > 0 & end(ChrGenes)+1000 < ChrLengths[i]]
    	  #print(length(ChrGenes)-length(TrimmedGenes[[i]]))
	}
  }
  Genes <-  unlist(GRangesList(Genes))
  return(Genes)
}



#Arguments <- c("/lustre/mib-cri/carrol09/Work/20120810_LewisS_AM_MMhmeDIP/bamFiles/EightBound.bwa.Realignedmm9_Processed.bam",getwd(),"/lustre/mib-cri/carrol09/Annotation/mm9/mm9Genes.bed","/lustre/mib-cri/carrol09/MyPipe/bedFiles/mm9.txt")
#ChromosomeLengths <-  "/lustre/mib-cri/carrol09/MyPipe/bedFiles/mm9.txt"
fileForCov <- Arguments[1]
DirectoryOut <- Arguments[2]
ChromosomeLengths <- Arguments[4]


ChrLims <- read.delim(ChromosomeLengths,sep="\t",comment.char="#",header=F)


library(GenomicRanges)
print(Arguments[3])

  if(any(grepl(".RData",Arguments[3]))){
    Genes <- get(load(Arguments[3]))
  }
  if(any(grepl(".bed",Arguments[3]))){
    GenesTemp <- read.delim(Arguments[3],sep="\t",comment.char="#",header=T)
    Genes <- Bed2GeneGRanges(GenesTemp,header=T,strand=T)
  }

elementMetadata(Genes)$Feature <- "Gene"
TrimmedGenes <- OnlyGenesInLims(Genes,ChrLims)


CoverageFeature(fileForCov,4000,1000,"TSS_AvCov_",TrimmedGenes,DirectoryOut)

