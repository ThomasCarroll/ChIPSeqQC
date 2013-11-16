##Get Arguments
Args <- commandArgs(trailingOnly = TRUE)

## Get the Working directory from the supplied argument
WkgDir <- getwd()
#WkgDir <- Args[2]
#getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
JobString <- "toysd2"
#JobString <- Args[3]

WkgDir <- getwd()

GetPipelinebase <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  
  PipelineBase <- ConfigFile[ConfigFile[,2] %in% "BaseLocation",3]
  return(PipelineBase)
}  
PipelineBase <- GetPipelinebase()

source(file.path(PipelineBase,"/RScripts/Workflow_Functions3.r"))



PipeLineLocations <- GetImportantLocations(WkgDir,"Config")


InFile <- Args[1]
OutFile <- Args[2]
RankColumn <- as.numeric(Args[3])
RankOrder <- Args[4]
NPeaks <- as.numeric(Args[5])



#InFile <- "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/Trial/20121114_MontoyaAR_DN_Hes6ChIP/Peaks/Macs/MCF7.bwa_Processed_peaks.bed"
#OutFile <- "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/Trial/20121114_MontoyaAR_DN_Hes6ChIP/Peaks/Macs/MCF7.bwa_Processed_peaksTop500.bed"
#RankColumn <- as.numeric(5)
#RankOrder <- "Rank"
#NPeaks <- as.numeric(500)


PeaksForMotifs <- TopPeaksByRank(InFile,OutFile,RankColumn,RankOrder,NPeaks)


