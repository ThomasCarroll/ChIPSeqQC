###This subscript deals with the samplesheet and finding locations.

##Get Arguments
Args <- commandArgs(trailingOnly = TRUE)

## Get the Working directory from the supplied argument
WkgDir <- getwd()
#WkgDir <- Args[2]
#getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
#JobString <- "toysd2"
#JobString <- Args[3]

JobString <- Args[1]


GetPipelinebase <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  
  PipelineBase <- ConfigFile[ConfigFile[,2] %in% "BaseLocation",3]
  return(PipelineBase)
}  
PipelineBase <- GetPipelinebase()


source(file.path(PipelineBase,"/RScripts/Workflow_Functions3.r"))
## Parse from config important locations

PipeLineLocations <- GetImportantLocations(WkgDir,"Config")
## Create the directory structure based on config
#CreateDirStruct(WkgDir,"Config",PipeLineLocations)
## Read SampleSheet if it exists..

RunConfigureAnnotationPipeline2(WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipeLineLocations,Config="Config")

PipeLineLocations <- GetImportantLocations(WkgDir,"Config")

if(!file.exists(file.path(PipeLineLocations@LocationsDir,"Lims_SampleLocations.txt"))){
    file.create(file.path(PipeLineLocations@LocationsDir,"Lims_SampleLocations.txt"))
}


SampleSheet <- GetSampleSheetOrDummy(WkgDir,"SampleSheet.csv")

## Parse what is local (FQs in Bam or FQ directory) or in Lims
SampleSheet <- GetLocal(SampleSheet,PipeLineLocations)
#LimsLog <- GetSamplesFromLims(SampleSheet)

# write.table(SampleSheet,"SetUp.csv",sep=",",row.names=F)
SampleSheet <- RefreshSampleSheet(SampleSheet,PipeLineLocations)




##Check for samples to be demultiplexed
SampleSheet <- CheckFromDemultiplex(PLs=PipeLineLocations,SS=SampleSheet)


#WriteLog(WhatsLocal,LimsLog,ToDemultiplex,PipeLineLocations)

SampleSheet <- RunSSfetchPipeline(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=300,PLs=PipeLineLocations,Config="Config")
### i need to add number of trys!!

BamsToRealign <- RunCheckGenomePipeline(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=300,PLs=PipeLineLocations,Config="Config")

SampleSheet <- RunSSfqfetchPipeline(SampleToGrab=BamsToRealign,SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=400,PLs=PipeLineLocations,Config="Config")

SampleSheet <- RunSSdeMultiplesPipeline(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=300,PLs=PipeLineLocations,Config="Config")

SampleSheet <- RunSSRealignmentPipeline(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=400,PLs=PipeLineLocations,Config="Config")

SampleSheet <- RunSSMergingPipeline(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipeLineLocations,Config="Config")

SampleSheet <- RunBamProcessPipeline(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipeLineLocations,Config="Config")

SampleSheet <- RunDownSamplePipeline(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipeLineLocations,Config="Config")


write.table(SampleSheet,"SampleSheet.csv",sep=",",row.names=F,quote=F)



write.table("Complete",file.path(PipeLineLocations@WorkFlowDir,paste(JobString,"_SampleSheetSetup.txt",sep="")),col.names=T,row.names=F,sep=",",quote=F)



