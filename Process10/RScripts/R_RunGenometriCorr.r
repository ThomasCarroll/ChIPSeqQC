
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
library(raster)
library(GenometriCorr)

S4Peaks <- read.delim(Args[1],comment.char="#",header=T)
PeakGranges1 <- Bed2GRanges(S4Peaks,header=T,strand=F)
S5Peaks <- read.delim(Args[2],comment.char="#",header=T)
PeakGranges2 <- Bed2GRanges(S5Peaks,header=T,strand=F)

ConfigFile <- readIniFile(file.path(getwd(),"Temp","config.ini"))
GENOME_PARAMETER <- ConfigFile[ConfigFile[,2] %in% "genome",3]

if(tolower(GENOME_PARAMETER) %in% tolower("GRCh37")){
	GENOME_LENGTHS <- "/lustre/mib-cri/carrol09/MyPipe/bedFiles/hg19.txt"
}
if(tolower(GENOME_PARAMETER) %in% tolower("HG18")){
	GENOME_LENGTHS <- "/lustre/mib-cri/carrol09/MyPipe/bedFiles/hg18.txt"
}
Genome_Lengths <- read.delim(GENOME_LENGTHS,sep="\t",header=F)
GenomeLengths <- as.numeric(as.vector(Genome_Lengths[,2]))
names(GenomeLengths) <- as.vector(Genome_Lengths[,1])

Peak1ToPeak2 <- GenometriCorrelation(PeakGranges1,PeakGranges2,awhole.only=T)
Peak1ToPeak2Res <- print(Peak1ToPeak2)
Res1 <- unlist(Peak1ToPeak2@.Data)
Res <- cbind(names(Res1),as.vector(Res1))
write.table(Res,Args[3],sep="\t",col.names=F,row.names=F)


