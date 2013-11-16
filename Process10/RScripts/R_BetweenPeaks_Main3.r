###This subscript deals with the samplesheet and processing Bamfiles and finding predicted fragment lengths.

##Get Arguments
Args <- commandArgs(trailingOnly = TRUE)

## Get the Working directory from the supplied argument
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
#
#WkgDir <- Args[2]
#JobString <- "toysd2"
JobString <- Args[1]

PeakCaller <- Args[3]

## Parse from config important locations
PipeLineLocations <- GetImportantLocations(WkgDir,"Config")

SampleSheet <- read.delim("SampleSheet.csv",sep=",",header=T)

SampleSheet <- RunBetweenPeaksPipeline(SampleSheet,WkgDir,JobString,MaxJobs=75,PipeLineLocations,"Config",Caller=PeakCaller)

#write.table(SampleSheet,"SampleSheet.csv",sep=",",row.names=F,quote=F)

write.table("Complete",file.path(PipeLineLocations@WorkFlowDir,paste(JobString,"_",PeakCaller,"_BetweenPeaks_Main_3.txt",sep="")),col.names=T,row.names=F,sep=",",quote=F)


