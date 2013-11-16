
## Now some of my own nonsense.

CreateDirStruct <- function(WkgDir=getwd(),ConfigDirectory="Config",PLs=PipeLineLocations){

  dir.create(PLs@BamDir,recursive=F,showWarnings=F)
  dir.create(PLs@FQDir,recursive=F,showWarnings=F)
  dir.create(PLs@FragLengthDir,recursive=F,showWarnings=F)
  dir.create(PLs@LocationsDir,recursive=F,showWarnings=F)
  dir.create(PLs@TempDir,recursive=F,showWarnings=F)
  dir.create(PLs@WorkFlowDir,recursive=F,showWarnings=F)
  
  system(paste("mkdir -p ",PLs@MacsDir,sep=""),wait=T)  
  system(paste("mkdir -p ",PLs@SicerDir,sep=""),wait=T)  
  system(paste("mkdir -p ",PLs@TPICsDir,sep=""),wait=T)  
  StripeOrNot <- GetStripingFromConfig(WkgDir,ConfigDirectory="Config")      
  if(StripeOrNot == "True"){
    system(paste("lfs setstripe ",PLs@BamDir,sep=""),wait=T)
    system(paste("lfs setstripe ",PLs@FQDir,sep=""),wait=T)
    system(paste("lfs setstripe ",PLs@MacsDir,sep=""),wait=T)  
    system(paste("lfs setstripe ",PLs@SicerDir,sep=""),wait=T) 
    system(paste("lfs setstripe ",PLs@TPICsDir,sep=""),wait=T)     
  } 
}

UnlockSampleSheet <- function(WkgDir){
  if(file.exists(file.path(WkgDir,"SampleSheet.LOCK"))){
    unlink(file.path(WkgDir,"SampleSheet.LOCK"))
  }

}


ReadAndLock <- function(ss,WkdDir,SAF=T,napTime=5){
  if(file.exists(gsub(".csv",".LOCK",ss))){
    while(file.exists(gsub(".csv",".LOCK",ss))){
      Sys.sleep(napTime)
    }
    write.table("Locked",gsub(".csv",".LOCK",ss))
    SampleSheet <- read.delim(ss,stringsAsFactors=SAF,sep=",")
  }else{
    write.table("Lcoked",gsub(".csv",".LOCK",ss))
    SampleSheet <- read.delim(ss,stringsAsFactors=SAF,sep=",")
  }
  return(SampleSheet)

}

WriteAndUnlock <- function(SampleSheet,ss){
   if(file.exists(gsub(".csv",".LOCK",ss))){
     write.table(SampleSheet,ss,sep=",",row.names=F,quote=F)
     unlink(gsub(".csv",".LOCK",ss))
   }else{
     write.table(SampleSheet,ss,sep=",",row.names=F,quote=F)
   }
}

GetUsernames <- function(){
  UserNameList <- vector("list",length=2)
  UserName = system("id -nu",wait=TRUE,intern=TRUE)
  MegaName <- paste("cri.camres.org\\\\",UserName,"@uk-cri-larc01",sep="")
  UserNameList[[1]] <-  UserName
  UserNameList[[2]] <-  MegaName
  return(UserNameList)
}

getRandString<-function(len=12){
  return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
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
  CovDirTemp <- ConfigFile[ConfigFile[,2] %in% "coveragedirectory",3]  
  setClass("ChIPDirLocations", representation(TempDir = "character",WkgDir = "character",LocationsDir = "character",BamDir = "character",FQDir = "character",WorkFlowDir="character",FragLengthDir="character",MacsDir="character",SicerDir="character",TPICsDir="character",CovDir="character"))
  PLLocations <- new("ChIPDirLocations",TempDir=TempDirTemp,WkgDir=WkgDirTemp,LocationsDir=LocationsDirTemp,BamDir=BamDirTemp,FQDir=FQDirTemp,WorkFlowDir=WorkFlowDirTemp,FragLengthDir=FragLengthDirTemp,MacsDir=MacsDirTemp,SicerDir=SicerDirTemp,TPICsDir=TPICsDirTemp,CovDir=CovDirTemp)
  return(PLLocations)
}


GetGenomeFromConfig <- function(WkgDir,ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  return(ConfigFile[ConfigFile[,"section"] == "Analysis Settings" & ConfigFile[,"name"] == "genome","value"])  
}

GetGenomeBuildFromConfig <- function(WkgDir,ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  genome <- ConfigFile[ConfigFile[,"section"] == "Analysis Settings" & ConfigFile[,"name"] == "genome","value"]
  genomeBuild <- ConfigFile[ConfigFile[,"section"] == "Genomes" & tolower(ConfigFile[,"name"]) == tolower(genome),"value"]
  return(genomeBuild)
}

GetChrLengthsFromConfig <- function(WkgDir,ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  genome <- ConfigFile[ConfigFile[,"section"] == "Analysis Settings" & ConfigFile[,"name"] == "genome","value"]
  ChrLengths <- ConfigFile[ConfigFile[,"section"] == "Chromosome Lengths" & tolower(ConfigFile[,"name"]) == tolower(genome),"value"]
  return(ChrLengths)
}

GetSequenceDictionaryFromConfig <- function(WkgDir,ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  genome <- ConfigFile[ConfigFile[,"section"] == "Analysis Settings" & ConfigFile[,"name"] == "genome","value"]
  SeqDict <- ConfigFile[ConfigFile[,"section"] == "Sequence Dictionary" & tolower(ConfigFile[,"name"]) == tolower(genome),"value"]
  return(SeqDict)
}

GetExcludedListFromConfig <- function(WkgDir,ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  genome <- ConfigFile[ConfigFile[,"section"] == "Analysis Settings" & ConfigFile[,"name"] == "genome","value"]
  ExclRegions <- ConfigFile[ConfigFile[,"section"] == "Excluded Regions" & tolower(ConfigFile[,"name"]) == tolower(genome),"value"]
  return(ExclRegions)
}

GetStripingFromConfig <- function(WkgDir,ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  striping <- ConfigFile[ConfigFile[,"section"] == "Analysis Settings" & ConfigFile[,"name"] == "stripe","value"]
  return(striping)
}


GetGenePosFromConfig <- function(WkgDir,ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  genome <- ConfigFile[ConfigFile[,"section"] == "Analysis Settings" & ConfigFile[,"name"] == "genome","value"]
  GenePos <- ConfigFile[ConfigFile[,"section"] == "Gene Positions" & tolower(ConfigFile[,"name"]) == tolower(genome),"value"]
  return(GenePos)
}

GetGeneSetsFromConfig <- function(WkgDir,ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  genome <- ConfigFile[ConfigFile[,"section"] == "Analysis Settings" & ConfigFile[,"name"] == "genome","value"]
  GeneSets <- ConfigFile[ConfigFile[,"section"] == "GeneSets" & tolower(ConfigFile[,"name"]) == tolower(genome),"value"]
  return(GeneSets)
}


GetQualityEncodingFromConfig <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  QE <- ConfigFile[ConfigFile[,"section"] == "Analysis Settings" & ConfigFile[,"name"] == "qualityencoding","value"]
  return(QE)
}


getMacsGenome <- function(WkgDir,ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  genome <- ConfigFile[ConfigFile[,"section"] == "Analysis Settings" & ConfigFile[,"name"] == "genome","value"]
  genomeformacs <- ConfigFile[ConfigFile[,"section"] == "Macs Parameters" & tolower(ConfigFile[,"name"]) == tolower(genome),"value"]
  return(genomeformacs)
}

getMacsmfold <- function(WkgDir,ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  mfold <- ConfigFile[ConfigFile[,"section"] == "Macs Parameters" & ConfigFile[,"name"] == "mfold","value"]
  return(mfold)
}

getMacsShiftSize <- function(WkgDir,ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  ShiftSize <- ConfigFile[ConfigFile[,"section"] == "Macs Parameters" & ConfigFile[,"name"] == "shiftsizedefault","value"]
  return(ShiftSize)
}

getSicerWindowSize <- function(WkgDir,ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  WindowSize <- ConfigFile[ConfigFile[,"section"] == "Sicer Parameters" & ConfigFile[,"name"] == "window","value"]
  return(WindowSize)
}

getSicerGapSize <- function(WkgDir,ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  GapSize <- ConfigFile[ConfigFile[,"section"] == "Sicer Parameters" & ConfigFile[,"name"] == "gapsize","value"]
  return(GapSize)
}

getFastxTrimmerExec  <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  exec <- ConfigFile[ConfigFile[,"section"] == "Executables" & ConfigFile[,"name"] == "fastxtrimmer","value"]
  return(exec)
}

getReadLengthMax  <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  exec <- ConfigFile[ConfigFile[,"section"] == "Analysis Settings" & ConfigFile[,"name"] == "readlengthmax","value"]
  return(exec)
}




getSicerExec  <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  exec <- ConfigFile[ConfigFile[,"section"] == "Executables" & ConfigFile[,"name"] == "sicer","value"]
  return(exec)
}

getSicerCustomExec  <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  exec <- ConfigFile[ConfigFile[,"section"] == "Custom Scripts" & ConfigFile[,"name"] == "sicer_cri_script","value"]
  return(exec)
}

getTPICsCustomExec  <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  exec <- ConfigFile[ConfigFile[,"section"] == "Custom Scripts" & ConfigFile[,"name"] == "tpics_cri_script","value"]
  return(exec)
}
getTPICsZetaCustomExec  <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  exec <- ConfigFile[ConfigFile[,"section"] == "Custom Scripts" & ConfigFile[,"name"] == "tpicszeta_cri_script","value"]
  return(exec)
}
getTPICsCreateCoverageCustomExec  <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  exec <- ConfigFile[ConfigFile[,"section"] == "Custom Scripts" & ConfigFile[,"name"] == "tpicscreatecoverage_cri_script","value"]
  return(exec)
}


getTPICSminSize <- function(WkgDir,ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  minsize <- ConfigFile[ConfigFile[,"section"] == "TPICs Parameters" & ConfigFile[,"name"] == "minsize","value"]
  return(minsize)
}

getTPICSwideSize  <- function(WkgDir,ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  wideSize <- ConfigFile[ConfigFile[,"section"] == "TPICs Parameters" & ConfigFile[,"name"] == "widesizeregion","value"]
  return(wideSize)
}


CallMotifsCheck <- function(Caller,WkgDir=getwd(),Config="Config"){
   Flags <- GetFlags(WkgDir,Config)
   LogicalCall <- slot(Flags,paste("Call",Caller,"Motifs",sep="")) == "Yes"
   return(LogicalCall)
}

CallProfilesCheck <- function(Caller,WkgDir=getwd(),Config="Config"){
   Flags <- GetFlags(WkgDir,Config)
   LogicalCall <- slot(Flags,paste("Call",Caller,"PeakProfile",sep="")) == "Yes"
   return(LogicalCall)
}

CallPeaksCheck <- function(Caller,WkgDir=getwd(),Config="Config"){
   Flags <- GetFlags(WkgDir,Config)
   LogicalCall <- slot(Flags,paste("Call",Caller,"Peaks",sep="")) == "Yes"
   return(LogicalCall)
}

AutoMergeCheck <- function(WkgDir=getwd(),Config="Config"){
   Flags <- GetFlags(WkgDir,Config)
   LogicalCall <- slot(Flags,"AutoMerging") == "Yes"
   return(LogicalCall)
}

CallBetweenPeaksCheck <- function(Caller,WkgDir=getwd(),Config="Config"){
   Flags <- GetFlags(WkgDir,Config)
   LogicalCall <- slot(Flags,paste("Call",Caller,"BetweenPeaks",sep="")) == "Yes"
   return(LogicalCall)
}

GetDupFlag <- function(WkgDir=getwd(),Config="Config"){
   Flags <- GetFlags(WkgDir,Config)
   LogicalCall <- slot(Flags,"RmDup")
   return(LogicalCall)
}

GetExcludedFlag <- function(WkgDir=getwd(),Config="Config"){
   Flags <- GetFlags(WkgDir,Config)
   LogicalCall <- slot(Flags,"ExcRegFlag")
   return(LogicalCall)
}

GetMapQFlag <- function(WkgDir=getwd(),Config="Config"){
   Flags <- GetFlags(WkgDir,Config)
   LogicalCall <- slot(Flags,"MapqFlag")
   if(as.numeric(LogicalCall) > 0){
     LogicalCall <- "True"
   }else{LogicalCall <- "False"}
   return(LogicalCall)
}



GetFlags <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  
  MapqFlagTemp <- ConfigFile[ConfigFile[,2] %in% "mapqfilter",3]
  ExcRegFlagTemp <- ConfigFile[ConfigFile[,2] %in% "useexcludedregionfilter",3]
  RmDupTemp <- ConfigFile[ConfigFile[,2] %in% "removeduplicates",3]
  CallMacsPeaksTemp <- ConfigFile[ConfigFile[,2] %in% "callmacspeaks",3]
  CallMacsMotifsTemp <- ConfigFile[ConfigFile[,2] %in% "callmacsmotifs",3]
  CallMacsPeakProfileTemp <- ConfigFile[ConfigFile[,2] %in% "callmacspeakprofile",3]  
  CallMacsBetweenPeaksTemp <- ConfigFile[ConfigFile[,2] %in% "callmacsbetweenpeaks",3] 
  CallSicerPeaksTemp <- ConfigFile[ConfigFile[,2] %in% "callsicerpeaks",3]
  CallSicerMotifsTemp <- ConfigFile[ConfigFile[,2] %in% "callsicermotifs",3]
  CallSicerPeakProfileTemp <- ConfigFile[ConfigFile[,2] %in% "callsicerpeakprofile",3]  
  CallTpicsPeaksTemp <- ConfigFile[ConfigFile[,2] %in% "calltpicspeaks",3]
  CallTpicsMotifsTemp <- ConfigFile[ConfigFile[,2] %in% "calltpicsmotifs",3]
  CallTpicsPeakProfileTemp <- ConfigFile[ConfigFile[,2] %in% "calltpicspeakprofile",3] 
  AutoMergingTemp <- ConfigFile[ConfigFile[,2] %in% "automerging",3]
  CallMacsBetweenPeaksTemp <- ConfigFile[ConfigFile[,2] %in% "callmacsbetweenpeaks",3] 
  CallSicerBetweenPeaksTemp <- ConfigFile[ConfigFile[,2] %in% "callsicerbetweenpeaks",3] 
  CallTPICsBetweenPeaksTemp <- ConfigFile[ConfigFile[,2] %in% "calltpicsbetweenpeaks",3]       
  
     
  setClass("ChIPFlags", representation(
  MapqFlag = "character",ExcRegFlag = "character",RmDup = "character",CallMacsPeaks = "character",CallMacsMotifs = "character",
  CallMacsPeakProfile="character",CallSicerPeaks="character",CallSicerMotifs="character",CallSicerPeakProfile="character",
  CallTPICsPeaks="character",CallTPICsMotifs="character",CallTPICsPeakProfile="character",AutoMerging="character",
  CallMacsBetweenPeaks="character",CallSicerBetweenPeaks="character",CallTPICsBetweenPeaks="character")
  )
  Flags <- new("ChIPFlags",
  MapqFlag=MapqFlagTemp,ExcRegFlag=ExcRegFlagTemp,RmDup=RmDupTemp,CallMacsPeaks=CallMacsPeaksTemp,CallMacsMotifs=CallMacsMotifsTemp,
  CallMacsPeakProfile=CallMacsPeakProfileTemp,CallSicerPeaks=CallSicerPeaksTemp,CallSicerMotifs=CallSicerMotifsTemp,
  CallSicerPeakProfile=CallSicerPeakProfileTemp,CallTPICsPeaks=CallTpicsPeaksTemp,CallTPICsMotifs=CallTpicsMotifsTemp,
  CallTPICsPeakProfile=CallTpicsPeakProfileTemp,AutoMerging=AutoMergingTemp,CallMacsBetweenPeaks=CallMacsBetweenPeaksTemp,
  CallSicerBetweenPeaks=CallSicerBetweenPeaksTemp,CallTPICsBetweenPeaks=CallTPICsBetweenPeaksTemp
  )
  return(Flags)
}

GetExecConfig <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  
  bwa <- ConfigFile[ConfigFile[,2] %in% "bwa",3]
  python <- ConfigFile[ConfigFile[,2] %in% "python",3]
  samtools <- ConfigFile[ConfigFile[,2] %in% "samtools",3]
  picard <- ConfigFile[ConfigFile[,2] %in% "picard",3]
  perl <- ConfigFile[ConfigFile[,2] %in% "perl",3]
  rsync <- ConfigFile[ConfigFile[,2] %in% "rsync",3]
  bedtools <- ConfigFile[ConfigFile[,2] %in% "bedtools",3]
  java <- ConfigFile[ConfigFile[,2] %in% "java",3]              
  rexec <- ConfigFile[ConfigFile[,2] %in% "rexec",3]  
  bigwig <- ConfigFile[ConfigFile[,2] %in% "bigwig",3]  
  macs <- ConfigFile[ConfigFile[,2] %in% "macs",3]  
  meme <- ConfigFile[ConfigFile[,2] %in% "meme",3]  
  ame <- ConfigFile[ConfigFile[,2] %in% "ame",3]  
  sicer <- ConfigFile[ConfigFile[,2] %in% "sicer",3]  
  tpicszeta <- ConfigFile[ConfigFile[,2] %in% "tpicszeta",3]   
  tpicscreatecoverage <- ConfigFile[ConfigFile[,2] %in% "tpicscreatecoverage",3]
  tpics <- ConfigFile[ConfigFile[,2] %in% "tpics",3]  
  gtftobed <- ConfigFile[ConfigFile[,2] %in% "gtftobed",3]
  
  setClass("ExecConfig", representation(
  bwa = "character",python = "character",samtools = "character",picard= "character",perl= "character",
  rsync = "character",bedtools = "character",java = "character",rexec = "character",bigwig="character",macs="character",ame="character",meme="character",
  sicer = "character",tpicszeta = "character",tpicscreatecoverage = "character",tpics = "character",gtftobed = "character"
  ))
  PLExec <- new("ExecConfig",
  bwa = bwa,python = python,samtools = samtools,picard= picard,perl= perl,
  rsync = rsync,bedtools = bedtools,java = java,rexec = rexec,bigwig=bigwig,macs=macs,ame=ame,meme=meme,
  sicer = sicer,tpicszeta = tpicszeta,tpicscreatecoverage = tpicscreatecoverage,tpics = tpics,gtftobed = gtftobed  
  )
  return(PLExec)
}

GetLibraryConfig <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  
  rlibs <- ConfigFile[ConfigFile[,2] %in% "rlibs",3]
  pythonlibs <- ConfigFile[ConfigFile[,2] %in% "pythonlibs",3]
  perllibs <- ConfigFile[ConfigFile[,2] %in% "perllibs",3]
  javalibs <- ConfigFile[ConfigFile[,2] %in% "javalibs",3]

  setClass("LibraryConfig", representation(
  rlibs = "character",pythonlibs = "character",perllibs = "character",javalibs = "character"
  ))
  PLLibrary <- new("LibraryConfig",
  rlibs = rlibs,pythonlibs = pythonlibs,perllibs = perllibs,javalibs = javalibs)
  return(PLLibrary)
}

getLibraryPath <- function(LibraryToRun,WkgDir=getwd(),Config="Config"){
   Libraries <- GetLibraryConfig(WkgDir,Config)
   LibraryPath <- slot(Libraries,paste(LibraryToRun,sep=""))
   return(LibraryPath)
}


getExecPath <- function(ExecToRun,WkgDir=getwd(),Config="Config"){
   Execs <- GetExecConfig(WkgDir,Config)
   ExecPath <- slot(Execs,paste(ExecToRun,sep=""))
   return(ExecPath)
}

GetTFDB <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  
  tfdb <- ConfigFile[ConfigFile[,2] %in% "tfdb",3]
  return(tfdb)
}  


DownSampleCheck <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  
  downsampleflag <- ConfigFile[ConfigFile[,2] %in% "downsample",3]
  return(downsampleflag=="Yes")
}  


GetPipelinesConfig <- function(WkgDir=getwd(),ConfigDirectory="Config"){

  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  
  mainpipeline <- ConfigFile[ConfigFile[,2] %in% "mainpipeline",3]
  bamfetchpipeline <- ConfigFile[ConfigFile[,2] %in% "bamfetchpipeline",3]
  checkgenomepipeline <- ConfigFile[ConfigFile[,2] %in% "checkgenomepipeline",3]
  fqfetchpipeline <- ConfigFile[ConfigFile[,2] %in% "fqfetchpipeline",3]
  alignpipeline <- ConfigFile[ConfigFile[,2] %in% "alignpipeline",3]
  mergingpipeline <- ConfigFile[ConfigFile[,2] %in% "mergingpipeline",3]
  bamprocesspipeline <- ConfigFile[ConfigFile[,2] %in% "bamprocesspipeline",3]
  bamprofilepipeline <- ConfigFile[ConfigFile[,2] %in% "bamprofilepipeline",3]
  macspeakcallpipeline <- ConfigFile[ConfigFile[,2] %in% "macspeakcallpipeline",3]
  peakprofilepipeline <- ConfigFile[ConfigFile[,2] %in% "peakprofilepipeline",3]
  motifpipeline <- ConfigFile[ConfigFile[,2] %in% "motifpipeline",3]
  betweenpeakspipeline  <- ConfigFile[ConfigFile[,2] %in% "betweenpeakspipeline",3]
  acrosspeakspipeline  <- ConfigFile[ConfigFile[,2] %in% "acrosspeakspipeline",3]
  sicerpeakcallpipeline  <- ConfigFile[ConfigFile[,2] %in% "sicerpeakcallpipeline",3]  
  tpicspeakcallpipeline  <- ConfigFile[ConfigFile[,2] %in% "tpicspeakcallpipeline",3]
  mainreportpipeline  <- ConfigFile[ConfigFile[,2] %in% "mainreportpipeline",3]
  downsamplepipeline  <- ConfigFile[ConfigFile[,2] %in% "downsamplepipeline",3]  
  configureannotationpipeline  <- ConfigFile[ConfigFile[,2] %in% "configureannotationpipeline",3]  
  
    
  setClass("PipelinesConfig", representation(
  mainpipeline = "character",bamfetchpipeline = "character",checkgenomepipeline = "character",fqfetchpipeline = "character",alignpipeline = "character",mergingpipeline = "character",bamprocesspipeline="character",bamprofilepipeline="character",macspeakcallpipeline="character",peakprofilepipeline="character",
  motifpipeline="character",betweenpeakspipeline = "character",acrosspeakspipeline = "character",tpicspeakcallpipeline = "character",sicerpeakcallpipeline = "character",mainreportpipeline = "character",downsamplepipeline = "character",configureannotationpipeline = "character"

  ))
  PLPipelines <- new("PipelinesConfig",
  mainpipeline = mainpipeline,bamfetchpipeline = bamfetchpipeline,checkgenomepipeline = checkgenomepipeline,fqfetchpipeline = fqfetchpipeline,alignpipeline = alignpipeline,mergingpipeline = mergingpipeline,bamprocesspipeline=bamprocesspipeline,
  bamprofilepipeline=bamprofilepipeline,macspeakcallpipeline=macspeakcallpipeline,peakprofilepipeline=peakprofilepipeline,motifpipeline=motifpipeline,betweenpeakspipeline = betweenpeakspipeline,acrosspeakspipeline = acrosspeakspipeline,tpicspeakcallpipeline = tpicspeakcallpipeline,sicerpeakcallpipeline = sicerpeakcallpipeline, mainreportpipeline = mainreportpipeline,downsamplepipeline = downsamplepipeline,configureannotationpipeline = configureannotationpipeline
  )
  return(PLPipelines)
}

getPipelinesPath <- function(PipelineToRun,WkgDir=getwd(),Config="Config"){
   Pipelines <- GetPipelinesConfig(WkgDir,Config)
   PipelinePath <- slot(Pipelines,paste(PipelineToRun,sep=""))
   return(PipelinePath)
}


GetSampleSheetOrDummy <- function(WkgDir=getwd(),SampleSheet="SampleSheet.csv"){
  if(file.exists(file.path(WkgDir,SampleSheet))){
    SampleSheet <- read.delim("SampleSheet.csv",sep=",",header=T,stringsAsFactors=F)
    if(length(grep("^bamFileName$",colnames(SampleSheet))) > 0){
        colnames(SampleSheet)[grep("^bamFileName$",colnames(SampleSheet))] <- "Source_File"
    }
    if(!length(grep("^BamLocation$",colnames(SampleSheet))) > 0){
              colnames(SampleSheet)[grep("^Location$",colnames(SampleSheet))] <- "BamLocation"
    }
    if(!length(grep("^FQLocation$",colnames(SampleSheet))) > 0){
              SampleSheet[,"FQLocation"] <- NA
    }
    if(!length(grep("^BarcodesFile$",colnames(SampleSheet))) > 0){
              SampleSheet[,"BarcodesFile"] <- NA
    }
    if(!length(grep("^Analysis_State$",colnames(SampleSheet))) > 0){
              SampleSheet[,"Analysis_State"] <- "Local_SetToRun"
    }


  }else{
    SampleSheet <- matrix(nrow=1,ncol=30)
    colnames(SampleSheet) <- c("GenomicsID","SampleName","Analysis_State","Run","Lane","BamLocation","FQLocation","BarcodesFile","Tissue","Factor","Antibody","Condition_1","Condition_2","Replicate","InputToUse","toMerge","Source_File","Processed_bamFileName","Original","delRand","Excluded","Filtered","Unique","DuplicationRate","InputFileNameToUse","derivedFrom","Macs_Fragment_Length","Sissr_Fragment_Length","Correlation_Fragment_Length","Coverage_Fragment_Length")
  }
  return(SampleSheet)
}

GetLocal <- function(SS=SampleSheet,PLs=PipeLineLocations){
#  SamplesWeHave <- NULL
#  SamplesFoundFromPrevious <- vector("character")
#  SamplesNotFoundPrevious <- vector("character")
#  SamplesInLocations <- SS[,"Source_File"]
  AsBams <-  dir(file.path(PLs@BamDir),pattern="*.bam")
  FQDirs <-   list.dirs(file.path(PLs@FQDir),recursive=F)
  AllFQs <- dir(FQDirs,full.names=T)
  for(i in 1:nrow(SS)){
    if(SS[i,"Source_File"] %in% AsBams | SS[i,"Source_File"] %in% AllFQs){
      print("Found one!")
    }else{
       SS[i,"Source_File"] <- "Local_Sample_Not_Found"
    }
    
  }
  
  if(any(SS[,"Source_File"] %in% "Local_Sample_Not_Found")){
    write.table(SS[SS[,"Source_File"] %in% "Local_Sample_Not_Found","SampleName"],file.path(PLs@LocationsDir,"Missing_Local_Samples.txt"),sep="\t",quote=F,row.names=F)
    system(paste("cat",file.path(PLs@LocationsDir,"Missing_Local_Samples.txt"),file.path(PLs@LocationsDir,"SamplenamesFromLims.txt"),sep=" "),wait=T)
    #if(file.info(file.path(PLs@LocationsDir,"SamplenamesFromLims.txt"))$size > 1){
    #  SamplesInLims <- as.vector(read.delim(file.path(PLs@LocationsDir,"SamplenamesFromLims.txt"),header=F)[,1])
    #}else{
    #  SamplesInLims <- ""
    #}
    #SamplesInLims <- as.vector(na.omit(unique(c(SS[SS[,"Source_File"] %in% "Local_Sample_Not_Found","SampleName"],SamplesInLims))))
    #write.table(SamplesInLims,file.path(PLs@LocationsDir,"SamplenamesFromLims.txt"),quote=F,row.names=F,col.names=F)
  }
  return(SS)     
}

#PythonScriptForLims=Scripts@LimsCall

GetSamplesFromLims <- function(SampleSheet=SampleSheet,PLs=PipeLineLocations,ProjectsFile="ProjectsFromLims.txt",SLXIDFile="SLXIDsFromLims.txt",SamplenamesFile="SamplenamesFromLims.txt",PythonForLims="/lustre/mib-cri/carrol09/python/PythonInstall/bin/python2.7",PythonScriptForLims="/lustre/mib-cri/carrol09/Work/MyPipe/Process10/PythonScripts/AnotherProcess.py"){
  CallLimsCommand <- paste(PythonForLims,PythonScriptForLims,PipeLineLocations@LocationsDir,ProjectsFile,SLXIDFile,SamplenamesFile,sep=" ")
  LimsLog <- system(CallLimsCommand,wait=T,intern=T)  
  ResFromLimsCall <- paste("From Lims", paste(LimsLog[grep("Found",LimsLog)],collapse="\n"),sep="\n")
  return(ResFromLimsCall)
}

GetBarcodes <- function(PLs=PipeLineLocations){
    BarcodesDirectory <- file.path(PipeLineLocations@LocationsDir,"barcodes")
    BarcodesFiles <-  dir(BarcodesDirectory,full.names=T)
    BarcodesFileNames <-  dir(BarcodesDirectory,full.names=F)
    Barcodes <- vector("list",length=length(BarcodesFileNames))
    if(length(Barcodes) != 0){
      for(i in 1:length(BarcodesFiles)){
       TempCodes <- read.delim(BarcodesFiles[i],sep="\t",header=F)
       Barcodes[[i]] <- TempCodes
      }
    }
    names(Barcodes) <- gsub(".txt","",BarcodesFileNames)
    return(Barcodes)
}

CheckFromDemultiplex <- function(PLs=PipeLineLocations,SS=SampleSheet){
  Barcodes <- GetBarcodes(PLs)
  TempBarcodes <- matrix(nrow=1,ncol=ncol(SampleSheet))
  colnames(TempBarcodes) <- colnames(SS)
  NotKnownSS <- SS[(SS[,"FQLocation"]) %in% "Location_Not_Known",]
  SS <-  SS[!(SS[,"FQLocation"]) %in% "Location_Not_Known",]
  if(nrow(SS) > 0){
  for(i in 1:nrow(SS)){
      ID <- SS[i,"GenomicsID"]
      if(length(grep("^SLX",ID)) > 0){
            ID <- gsub("\\..*","",ID)
      }
      if(any(names(Barcodes) %in% ID)){
          SS[i,"Analysis_State"] <- "Parent_Of_Sample"
          BarcodesForSamples <- Barcodes[names(Barcodes) %in% ID]          
          for(k in 1:nrow(BarcodesForSamples[[1]])){
            EvenTemperRow <- vector(length=ncol(SampleSheet)) 
            names(EvenTemperRow) <- colnames(SS)
            BarcodeInfo <- as.vector(BarcodesForSamples[[1]][k,])
            EvenTemperRow["GenomicsID"] <- paste(SS[i,"GenomicsID"],BarcodeInfo[,2],sep="")
            EvenTemperRow["Analysis_State"] <- "RunMe"            
            EvenTemperRow["Source_File"] <- SS[i,"GenomicsID"]
            EvenTemperRow["SampleName"] <- as.character(as.vector(BarcodeInfo[,2]))
            EvenTemperRow[c("Run","Lane")] <- SS[i,c("Run","Lane")]
            EvenTemperRow[c("FQLocation","BamLocation")] <- "Premultiplex"
            EvenTemperRow[EvenTemperRow == "FALSE"] <- NA            
            TempBarcodes <- rbind(TempBarcodes,EvenTemperRow)
          }
          SS[i,"BarcodesFiles"] <- dir(file.path(PLs@LocationsDir,"barcodes"),full.names=T)[dir(file.path(PLs@LocationsDir,"barcodes")) %in% paste(names(BarcodesForSamples),".txt",sep="")]
      }
  }
  }
  SS <- rbind(SS,NotKnownSS,TempBarcodes)
  SS <- SS[!is.na(SS[,1]),]
  return(SS)
}


WriteLog <- function(WhatsLocal=WhatsLocal,LimsLog=LimsLog,ToDemultiplex=ToDemultiplex,PLs=PipeLineLocations){
  cat(WhatsLocal,LimsLog,"\n\n",ToDemultiplex,"\n",file=file.path(PLs@LocationsDir,"SamplesLog.txt"),append = TRUE)
}

GetWorkFlowConfig <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  
  MetaVersion <- ConfigFile[ConfigFile[,2] %in% "metaversion",3]
  XSI <- ConfigFile[ConfigFile[,2] %in% "xsi",3]
  SchemaLocation <- ConfigFile[ConfigFile[,2] %in% "schemalocation",3]
  Mode <- ConfigFile[ConfigFile[,2] %in% "mode",3]
  Pipeline <- ConfigFile[ConfigFile[,2] %in% "pipeline",3]
  TaskDirectories <- ConfigFile[ConfigFile[,2] %in% "taskdirectories",3]
  TempDirectory <- ConfigFile[ConfigFile[,2] %in% "tempdirectory",3]
  SummaryFile <- ConfigFile[ConfigFile[,2] %in% "summaryfile",3]
  SummaryErrors <- ConfigFile[ConfigFile[,2] %in% "summaryerrors",3]
  Queue <- ConfigFile[ConfigFile[,2] %in% "queue",3]
  Executable <- ConfigFile[ConfigFile[,2] %in% "executable",3]
  setClass("WorkFlowConfig", representation(MetaVersion = "character",XSI = "character",SchemaLocation = "character",Mode= "character",Pipeline= "character",
  TaskDirectories = "character",TempDirectory = "character",SummaryFile = "character",SummaryErrors = "character",
  Queue = "character",Executable = "character"
  ))
  PLWorkFlow <- new("WorkFlowConfig",MetaVersion=MetaVersion,XSI=XSI,SchemaLocation=SchemaLocation,Mode=Mode,Pipeline=Pipeline,
  TaskDirectories = TaskDirectories,TempDirectory = TempDirectory,SummaryFile = SummaryFile,SummaryErrors = SummaryErrors,
  Queue = Queue,Executable = Executable  
  )
  return(PLWorkFlow)
}


getWorkflowParam <- function(WorkflowToRun,WkgDir=getwd(),Config="Config"){
   WorkflowParams <- GetWorkFlowConfig(WkgDir,Config)
   WorkflowParam <- slot(WorkflowParams,paste(WorkflowToRun,sep=""))
   return(WorkflowParam)
}

GetPipelinebase <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  
  PipelineBase <- ConfigFile[ConfigFile[,2] %in% "BaseLocation",3]
  return(PipelineBase)
}  

MakeWorkFlowXML<- function(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=75){
  require(XML)      
  WorkFlowConfig <- GetWorkFlowConfig(WkgDir=WkgDir,ConfigDirectory=Config)
  PLs <- GetImportantLocations(WkgDir,"Config")
  MetaVersion=WorkFlowConfig@MetaVersion                                           
  XSI=WorkFlowConfig@XSI
  SchemaLocation=WorkFlowConfig@SchemaLocation
  Mode = WorkFlowConfig@Mode
  TaskDirectories = unlist(strsplit(WorkFlowConfig@TaskDirectories,";"))
  TempDirectory =    WorkFlowConfig@TempDirectory
  SummaryErrors = WorkFlowConfig@SummaryErrors
  #lsfOutputDirectory = WorkFlowConfig@lfsOutputDirectory
  queue = WorkFlowConfig@Queue

  Pipeline = Pipeline
  maxJobs = maxJobs
  maxCpuResources = 2
  MainVariables <- Variables
  names(MainVariables) <- names(Variables)

  dir.create(file.path(PLs@WorkFlowDir,paste(PipeName,JobString,sep="")),showWarnings = FALSE)
  SummaryFile   =    file.path(PLs@WorkFlowDir,paste(PipeName,JobString,sep=""),"SummaryFile.txt")
  lsfOutputDirectory = file.path(PLs@WorkFlowDir,paste(PipeName,JobString,sep=""))

MetaXML <- newXMLDoc()
namespaceVariable <- c(MetaVersion,XSI)
names(namespaceVariable) <- c("meta","xsi")
attrs <- SchemaLocation
names(attrs) <- "xsi:schemaLocation"
GlobalNode <- newXMLNode("meta:metadata",namespace=namespaceVariable,attrs=attrs,suppressNamespaceWarning=T,parent = MetaXML)
  PipelineNode <- newXMLNode("pipeline",Pipeline)
  tempDirectoryNode <- newXMLNode("tempDirectory",TempDirectory)
  errorsonly <- SummaryErrors
  names(errorsonly) <- "errorsOnly"
  SummaryFileNode <- newXMLNode("summaryFile",attrs=errorsonly,SummaryFile)
  ModeNode <- newXMLNode("mode",Mode)
  TaskNode <- newXMLNode("taskDirectories")
  for(i in 1:length(TaskDirectories)){
    newXMLNode("directory",TaskDirectories[i],parent=TaskNode)
  }
  ExecConfigNode <- newXMLNode("executionConfiguration")
  modeattrs <- "lsf"
  names(modeattrs) <- "mode"
    executionNode1 <- newXMLNode("execution",attrs=modeattrs,parent=ExecConfigNode)
      LSFOutputNode <- newXMLNode("lsfOutputDirectory",lsfOutputDirectory,parent=executionNode1)
      queueNode <- newXMLNode("queue",queue,parent=executionNode1)
      MaxJobsNode <- newXMLNode("maximumSubmittedJobs",Pipeline,parent=executionNode1)
      localattrs <- "local"
      names(localattrs) <- "mode"
    executionNode2 <- newXMLNode("execution",attrs=localattrs,parent=ExecConfigNode)
      MaxJobsNode <- newXMLNode("maxCpuResources",Pipeline,parent=executionNode2)
   MainVariableNode <- newXMLNode("variables")
   for(i in 1:length(MainVariables)){
    newXMLNode(names(MainVariables)[i],MainVariables[i],parent=MainVariableNode)
   }
   if(length(Specialisations) > 0){
   MainSpecialNode <- newXMLNode("specialisations")
   SpecialisationNodeList <- vector("list",length=length(Specialisations))
   for(i in 1:length(Specialisations)){
   identifierattrs <- c(names(Specialisations)[i],"true")
   names(identifierattrs) <- c("identifier","active")
      SpecialisationNodeList[[i]] <- newXMLNode("specialisation",attrs=identifierattrs)
      VectorNodeList <- vector("list",length=length(SpecialisationNodeList[[i]]))  
      TempVariableNode <- newXMLNode("variables",parent=SpecialisationNodeList[[i]])
         for(k in 1:length(Specialisations[[i]])){
            newXMLNode(names(Specialisations[[i]])[k],Specialisations[[i]][k],parent=TempVariableNode)
          }
   }
   MainSpecialNode <- addChildren(MainSpecialNode,kids=SpecialisationNodeList,cdata=F)
   
  MetafileMainNode <-  list(PipelineNode,tempDirectoryNode,SummaryFileNode,ModeNode,TaskNode,ExecConfigNode,MainVariableNode,MainSpecialNode)
  }else{
   MetafileMainNode <-  list(PipelineNode,tempDirectoryNode,SummaryFileNode,ModeNode,TaskNode,ExecConfigNode,MainVariableNode)  
  }
  GlobalNode2 <- addChildren(GlobalNode,kids=MetafileMainNode,cdata=F)

  return(MetaXML)
}

MakeWorkFlowXML2<- function(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=75){
  require(XML)

  WorkFlowConfig <- GetWorkFlowConfig(WkgDir=WkgDir,ConfigDirectory=Config)
  PLs <- GetImportantLocations(WkgDir,"Config")
  MetaVersion=WorkFlowConfig@MetaVersion                                           
  XSI=WorkFlowConfig@XSI
  SchemaLocation=WorkFlowConfig@SchemaLocation
  Mode = WorkFlowConfig@Mode
  TaskDirectories = unlist(strsplit(WorkFlowConfig@TaskDirectories,";"))

  system(paste("mkdir -p ",file.path(WorkFlowConfig@TempDirectory,PipeName),sep=""))      
  TempDirectory =    file.path(WorkFlowConfig@TempDirectory,PipeName)


  SummaryErrors = WorkFlowConfig@SummaryErrors
  #lsfOutputDirectory = WorkFlowConfig@lfsOutputDirectory
  queue = WorkFlowConfig@Queue

  Pipeline = Pipeline
  maxJobs = maxJobs
  maxCpuResources = 2
  MainVariables <- Variables
  names(MainVariables) <- names(Variables)

  dir.create(file.path(PLs@WorkFlowDir,paste(PipeName,JobString,sep="")),showWarnings = FALSE)
  SummaryFile   =    file.path(PLs@WorkFlowDir,paste(PipeName,JobString,sep=""),"SummaryFile.txt")
  lsfOutputDirectory = file.path(PLs@WorkFlowDir,paste(PipeName,JobString,sep=""))

MetaXML <- newXMLDoc()
namespaceVariable <- c(MetaVersion,XSI)
names(namespaceVariable) <- c("meta","xsi")
attrs <- SchemaLocation
names(attrs) <- "xsi:schemaLocation"
GlobalNode <- newXMLNode("meta:metadata",namespace=namespaceVariable,attrs=attrs,suppressNamespaceWarning=T,parent = MetaXML)
  PipelineNode <- newXMLNode("pipeline",Pipeline)
  tempDirectoryNode <- newXMLNode("tempDirectory",TempDirectory)
  errorsonly <- SummaryErrors
  names(errorsonly) <- "errorsOnly"
  SummaryFileNode <- newXMLNode("summaryFile",attrs=errorsonly,SummaryFile)
  ModeNode <- newXMLNode("mode",Mode)
  TaskNode <- newXMLNode("taskDirectories")
  for(i in 1:length(TaskDirectories)){
    newXMLNode("directory",TaskDirectories[i],parent=TaskNode)
  }
  ExecConfigNode <- newXMLNode("executionConfiguration")
  modeattrs <- "lsf"
  names(modeattrs) <- "mode"
    executionNode1 <- newXMLNode("execution",attrs=modeattrs,parent=ExecConfigNode)
      LSFOutputNode <- newXMLNode("lsfOutputDirectory",lsfOutputDirectory,parent=executionNode1)
      queueNode <- newXMLNode("queue",queue,parent=executionNode1)
      MaxJobsNode <- newXMLNode("maximumSubmittedJobs",maxJobs,parent=executionNode1)
      localattrs <- "local"
      names(localattrs) <- "mode"
    executionNode2 <- newXMLNode("execution",attrs=localattrs,parent=ExecConfigNode)
      MaxJobsNode <- newXMLNode("maxCpuResources",maxCpuResources,parent=executionNode2)
   MainVariableNode <- newXMLNode("variables")
   for(i in 1:length(MainVariables)){
    newXMLNode(names(MainVariables)[i],MainVariables[i],parent=MainVariableNode)
   }
   if(length(Specialisations) > 0){
   MainSpecialNode <- newXMLNode("specialisations")
   SpecialisationNodeList <- vector("list",length=length(Specialisations))
   for(i in 1:length(Specialisations)){
   identifierattrs <- c(names(Specialisations)[i],"true")
   names(identifierattrs) <- c("identifier","active")
      SpecialisationNodeList[[i]] <- newXMLNode("specialisation",attrs=identifierattrs)
      VectorNodeList <- vector("list",length=length(SpecialisationNodeList[[i]]))  
      TempVariableNode <- newXMLNode("variables",parent=SpecialisationNodeList[[i]])
         for(k in 1:length(Specialisations[[i]])){
            newXMLNode(names(Specialisations[[i]])[k],Specialisations[[i]][k],parent=TempVariableNode)
          }
   }
   MainSpecialNode <- addChildren(MainSpecialNode,kids=SpecialisationNodeList,cdata=F)
   
  MetafileMainNode <-  list(PipelineNode,tempDirectoryNode,SummaryFileNode,ModeNode,TaskNode,ExecConfigNode,MainVariableNode,MainSpecialNode)
  }else{
   MetafileMainNode <-  list(PipelineNode,tempDirectoryNode,SummaryFileNode,ModeNode,TaskNode,ExecConfigNode,MainVariableNode)  
  }
  GlobalNode2 <- addChildren(GlobalNode,kids=MetafileMainNode,cdata=F)

  return(MetaXML)
}




RefreshSampleSheet <- function(SS=SampleSheet,PLs=PipeLineLocations){
  SS_Samples <- unique(as.character(SS[,1]))
#  SS_Samples_SLXIDS <-  gsub("\\..*","",SS_Samples[grep("^SLX",SS_Samples)])
 # names(SS_Samples_SLXIDS) <- SS_Samples[grep("^SLX",SS_Samples)]
  #SS_Samples_NoSLXIDS <-  SS_Samples[!grep("^SLX",SS_Samples)]
 # names(SS_Samples_NoSLXIDS) <- SS_Samples[!grep("^SLX",SS_Samples)]

#  LocalSamples <- read.delim()
  if(file.info(file.path(PLs@LocationsDir,"Lims_SampleLocations.txt"))$size > 0){
  LimsSamples <- read.delim(file.path(PLs@LocationsDir,"Lims_SampleLocations.txt"),sep="\t",header=F)
  }else{
  LimsSamples <- matrix(nrow=2,ncol=9)
  }
  
  NewLimsFrame <- cbind(paste(LimsSamples[,1],".",LimsSamples[,4],".s_",LimsSamples[,3],sep=""),as.vector(LimsSamples[,1]),as.vector(LimsSamples[,2]),paste("Run-",LimsSamples[,4],sep=""),paste("Lane-",LimsSamples[,3],sep=""),file.path(LimsSamples[,5],LimsSamples[,6]),file.path(LimsSamples[,7],LimsSamples[,8]),as.vector(LimsSamples[,9]))
  UpdateLimsFrame <- NewLimsFrame[NewLimsFrame[,1] %in% SS_Samples,]
  NewLimsFrame <- NewLimsFrame[!NewLimsFrame[,1] %in% SS_Samples,,drop=F]
  NewLimsFrame <- NewLimsFrame[match(unique(NewLimsFrame[,1]),NewLimsFrame[,1]),,drop=F]
  
  ##Make NewSS
  GenomicsID <- c(SS[match(SS[,1],SS_Samples),"GenomicsID"],as.vector(NewLimsFrame[,1]))
  SampleName <- c(SS[match(SS[,1],SS_Samples),"SampleName"],as.vector(NewLimsFrame[,3]))
  Run <- c(SS[match(SS[,1],SS_Samples),"Run"],as.vector(NewLimsFrame[,4]))
  Lane <- c(SS[match(SS[,1],SS_Samples),"Lane"],as.vector(NewLimsFrame[,5]))
  BamLocation <-  c(SS[match(SS[,1],SS_Samples),"BamLocation"],as.vector(NewLimsFrame[,6]))
  FQLocation <-  c(SS[match(SS[,1],SS_Samples),"FQLocation"],as.vector(NewLimsFrame[,7]))
  Analysis_State <- c(SS[match(SS[,1],SS_Samples),"Analysis_State"],as.vector(NewLimsFrame[,8]))
   
  BarcodesFilesTemp <- c(SS[match(SS[,1],SS_Samples),"BarcodesFile"])
  BarcodesFiles <- c(BarcodesFilesTemp,rep(NA,length(BamLocation)-length(BarcodesFilesTemp)))
#  Analysis_State <- c(Analysis_StateTemp,rep(NA,length(BamLocation)-length(Analysis_StateTemp))) 
  if(nrow(SS) > 1){   
    TempMetaData <- as.matrix(SS[,colnames(SS) %in% c("Tissue","Factor","Antibody","Condition_1","Condition_2","Replicate","InputToUse","toMerge")],ncol=8)
  }else{
    TempMetaData <- matrix(SS[,colnames(SS) %in% c("Tissue","Factor","Antibody","Condition_1","Condition_2","Replicate","InputToUse","toMerge")],ncol=8)
  }
  colnames(TempMetaData) <- c("Tissue","Factor","Antibody","Condition_1","Condition_2","Replicate","InputToUse","toMerge")
  Metadata <- matrix(nrow=length(GenomicsID),ncol=ncol(TempMetaData))  
  Metadata[1:nrow(SS),] <- as.matrix(TempMetaData)
  colnames(Metadata) <- colnames(SS)[colnames(SS) %in% c("Tissue","Factor","Antibody","Condition_1","Condition_2","Replicate","InputToUse","toMerge")]
  Source_File <- c(SS[match(SS[,1],SS_Samples),"Source_File"],rep("NA",length(as.vector(NewLimsFrame[,1]))))
  Processed_bamFileName <- c(SS[match(SS[,1],SS_Samples),"Processed_bamFileName"],rep("NA",length(as.vector(NewLimsFrame[,1]))))
  if(nrow(SS) > 1){
   TempResdata <- as.matrix(SS[,colnames(SS) %in% c("Original","delRand","Excluded","Filtered","Unique","DuplicationRate","InputFileNameToUse","derivedFrom","Macs_Fragment_Length","Sissr_Fragment_Length","Correlation_Fragment_Length","Coverage_Fragment_Length")],ncol=12)
  }else{  
   TempResdata <- matrix(SS[,colnames(SS) %in% c("Original","delRand","Excluded","Filtered","Unique","DuplicationRate","InputFileNameToUse","derivedFrom","Macs_Fragment_Length","Sissr_Fragment_Length","Correlation_Fragment_Length","Coverage_Fragment_Length")],ncol=12)  
  }
    colnames(TempResdata) <- c("Original","delRand","Excluded","Filtered","Unique","DuplicationRate","InputFileNameToUse","derivedFrom","Macs_Fragment_Length","Sissr_Fragment_Length","Correlation_Fragment_Length","Coverage_Fragment_Length")


  Resdata <- matrix(nrow=length(GenomicsID),ncol=ncol(TempResdata))
  Resdata[1:nrow(SS),] <- as.matrix(TempResdata)
  colnames(Resdata) <- colnames(SS)[colnames(SS) %in% c("Original","delRand","Excluded","Filtered","Unique","DuplicationRate","InputFileNameToUse","derivedFrom","Macs_Fragment_Length","Sissr_Fragment_Length","Correlation_Fragment_Length","Coverage_Fragment_Length")]
  
  CleanedSheet <- cbind(GenomicsID,SampleName,Run,Lane,Analysis_State,BamLocation,FQLocation,BarcodesFiles,Metadata,Source_File,Processed_bamFileName,Resdata)
  CleanedSheet <- CleanedSheet[!is.na(CleanedSheet[,"GenomicsID"]),]
  GoodNames <- GetUsernames()
  CleanedSheet[,"SampleName"] <- gsub(" ","_",CleanedSheet[,"SampleName"])
  CleanedSheet[,"BamLocation"] <- gsub("solexa/solexa","solexa",CleanedSheet[,"BamLocation"])
  CleanedSheet[,"FQLocation"] <- gsub("solexa/solexa","solexa",CleanedSheet[,"FQLocation"])
  CleanedSheet[is.na(CleanedSheet[,"FQLocation"]) | CleanedSheet[,"FQLocation"] %in% "/","FQLocation"] <-  "Location_Not_Known"
  CleanedSheet[is.na(CleanedSheet[,"BamLocation"]) | CleanedSheet[,"BamLocation"] %in% "/" | CleanedSheet[,"BamLocation"] %in% "NA/NA","BamLocation"] <-  "Location_Not_Known"
  
  CleanedSheet <- CleanedSheet[!apply(CleanedSheet,1,function(x)all(is.na(x))),]
  if(nrow(UpdateLimsFrame) > 0){
    for(i in 1:nrow(UpdateLimsFrame)){
        if(!is.na(UpdateLimsFrame[i,6]) & !UpdateLimsFrame[i,6] %in% "NA"){
            CleanedSheet[CleanedSheet[,1] %in% UpdateLimsFrame[i,1],"BamLocation"] <- gsub("solexa/solexa","solexa",UpdateLimsFrame[i,6])
        }  
        if(!is.na(UpdateLimsFrame[i,7]) & !UpdateLimsFrame[i,7] %in% "NA"){
            CleanedSheet[CleanedSheet[,1] %in% UpdateLimsFrame[i,1],"FQLocation"] <- gsub("solexa/solexa","solexa",UpdateLimsFrame[i,7])
        }            
          
        }
    }
    CleanedSheet <- CleanedSheet[!CleanedSheet[,1] == "NA.NA.s_NA",]
    CleanedSheet[,"BamLocation"] <- gsub("archive.cri.camres.org",GoodNames[[2]],as.vector(CleanedSheet[,"BamLocation"]),fixed=T)
    CleanedSheet[,"FQLocation"] <- gsub("archive.cri.camres.org",GoodNames[[2]],as.vector(CleanedSheet[,"FQLocation"]),fixed=T)
    CleanedSheet[,"BamLocation"] <- gsub("/sequencing/","/ArchiveDataFS1/Projects/bioinformatics-sequencing/",as.vector(CleanedSheet[,"BamLocation"]),fixed=T)
    CleanedSheet[,"FQLocation"] <- gsub("/sequencing/","/ArchiveDataFS1/Projects/bioinformatics-sequencing/",as.vector(CleanedSheet[,"FQLocation"]),fixed=T)

    for(i in 1:nrow(CleanedSheet)){
      AlternateName <- paste(unlist(strsplit(CleanedSheet[i,1],"\\."))[1],toupper(CleanedSheet[i,2]),unlist(strsplit(CleanedSheet[i,1],"\\."))[2],unlist(strsplit(CleanedSheet[i,1],"\\."))[3],sep=".")
      if(CleanedSheet[i,1] != gsub("\\.bwa\\.homo_sapiens\\.bam","",basename(CleanedSheet[i,"BamLocation"])) & AlternateName == gsub("\\.bwa\\.homo_sapiens\\.bam","",basename(CleanedSheet[i,"BamLocation"]))){
        CleanedSheet[i,1] <- AlternateName
      }
    }
    CleanedSheet <- CleanedSheet[match(unique(CleanedSheet[,1]),CleanedSheet[,1]),]



  return(CleanedSheet)
}


RunSSfetchPipeline <- function(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipeLineLocations,Config="Config"){
  
  Pipeline <- getPipelinesPath("bamfetchpipeline")
  workflowExec <- getWorkflowParam("Executable")
  javaExec <- getExecPath("java")
  mode <- tolower(getWorkflowParam("Mode"))
  PipelineBase <- GetPipelinebase()
  rsyncExec <- getExecPath("rsync")
  bamDir <- PLs@BamDir
  
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")
  
  
  Variables  <- c(file.path(WkgDir,""),rsyncExec,bamDir)
  names(Variables) <- c("WorkingDirectory","rsync","bamDir")
  #& SampleSheet[,"Analysis_State"] %in% "RunMe" 
  ToFetch <- SampleSheet[SampleSheet[,"BamLocation"] != "Location_Not_Known"  & SampleSheet[,"Analysis_State"] %in% "RunMe" & SampleSheet[,"BamLocation"] != "Premultiplex" & !is.na(SampleSheet[,"BamLocation"]) & SampleSheet[,"BamLocation"] != "/" & SampleSheet[,"BamLocation"] != "NA",c("GenomicsID","BamLocation")]
  if(length(ToFetch) > 0){
      Specialisations  <- vector("list",length=nrow(ToFetch))
      for(i in 1:length(Specialisations)){
        Specialisations[[i]] <- c(as.vector(ToFetch[i,"BamLocation"]),as.vector(ToFetch[i,"GenomicsID"]))
        names(Specialisations[[i]]) <- c("BamToGet","Bam_Name")  
      } 
      names(Specialisations) <-  paste(as.vector(ToFetch[,"GenomicsID"]),"_fetch",sep="") 
      #Pipeline = "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/BamGettingPipeline.xml"  
      PipeName <- "SS_Fetch"
      SSfetch_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=75)
      saveXML(SSfetch_WfMeta,file=file.path(PLs@WorkFlowDir,"SSfetch.xml"))
      system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,"SSfetch.xml"),sep=""),wait=T)
  }
  BamFilesReceived <- dir(path=PLs@BamDir,pattern="*.bam$")
  if(any(grep("Processed",BamFilesReceived))){BamFilesReceived <- BamFilesReceived[-grep("Processed",BamFilesReceived)]}
  if(any(grep("Realign",BamFilesReceived))){BamFilesReceived <- BamFilesReceived[-grep("Realign",BamFilesReceived)]}
  
  #BamID <- gsub(".bam","",gsub("\\.bwa.*","",BamFilesReceived))
  BamID <- BamFilesReceived
  for(i in 1:nrow(SampleSheet)){
  	if(!any(c(grep(".fq.gz$",SampleSheet[i,"Source_File"]),grep(".fq$",SampleSheet[i,"Source_File"]),grep("Processed",SampleSheet[i,"Source_File"]),grep("Realign",SampleSheet[i,"Source_File"])))){
      if(any(BamID %in% basename(SampleSheet[i,"BamLocation"]))){
  			SampleSheet[i,"Source_File"] <- BamID[BamID %in% basename(SampleSheet[i,"BamLocation"])]
  		}
  	}
  }

    SampleSheet <- SampleSheet[!apply(SampleSheet,1,function(x)all(is.na(x))),]

return(SampleSheet)
}



RunCheckGenomePipeline <- function(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipelineLocations,Config="Config"){

  Pipeline <- getPipelinesPath("checkgenomepipeline")
  workflowExec <- getWorkflowParam("Executable")
  javaExec <- getExecPath("java")
  mode <- tolower(getWorkflowParam("Mode"))
  pipelineBase <- GetPipelinebase()
  pythonExec <- getExecPath("python")
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")
  pythonlibs <- getLibraryPath("pythonlibs")
  bamDir <- PLs@BamDir
    

   
  genome <- GetGenomeFromConfig(WkgDir,ConfigDirectory=Config)
  Variables  <- c(file.path(WkgDir,""),genome,pythonExec,pipelineBase,pythonlibs,bamDir)
  names(Variables) <- c("WorkingDirectory","genome","python","pipelineBase","pythonlibs","bamDir")
  
  ToGenomify <- SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe" & SampleSheet[,"Source_File"] != "NA" & !is.na(SampleSheet[,"Source_File"]),"Source_File"]
  ToGenomify <- gsub(".bam","",ToGenomify[grep("*.bam",ToGenomify)])
  
  Specialisations  <- list()
  if(length(ToGenomify) > 0){
    for(i in 1:length(ToGenomify)){
      Specialisations[[i]] <- c(as.vector(ToGenomify[i]))
      names(Specialisations[[i]]) <- c("BamName")  
    } 
    names(Specialisations) <-  paste(as.vector(gsub("\\.bwa.*","",ToGenomify)),"_genomeCheck",sep="") 
  
#    Pipeline = "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/MultiGenomeGetterPipe.xml"  
    PipeName <- "SS_GenomeCheck"
    SSCheckGenome_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=75)
    saveXML(SSCheckGenome_WfMeta,file=file.path(PLs@WorkFlowDir,"SSCheckGenome.xml"))
    system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,"SSCheckGenome.xml"),sep=""),wait=T)
  }
    files <- dir(path=PLs@BamDir,pattern="*.info",full.names = T)
    Genomes <- vector("character",length=length(files))
    fileNames <- gsub(".info","",basename(files))
    SamplesFQToGrab <- vector()
    if(length(files) > 0){
    for(i in 1:length(files)){
  	   Genomes[i] <- as.vector(read.delim(files[i],header=F,sep="\t"))
    }
    Genomes <- as.vector(unlist(Genomes))
    GenomeMat <- cbind(paste(fileNames,".bam",sep=""),Genomes)
   
    GenomesToAlign <- GenomeMat[!tolower(GenomeMat[,2]) %in% tolower(genome),1]
    
    if(length(GenomesToAlign) > 0){
     SamplesFQToGrab <- vector("character",length=length(GenomesToAlign))
     for(i in 1:length(GenomesToAlign)){
      if(any(SampleSheet[,"Source_File"] %in% GenomesToAlign[i])){
  	   SamplesFQToGrab[i] <- as.vector(SampleSheet[SampleSheet[,"Source_File"] %in% GenomesToAlign[i],1])
       if(length(as.vector(SampleSheet[SampleSheet[,"Source_File"] %in% GenomesToAlign[i],1])) > 1){
       print(as.vector(SampleSheet[SampleSheet[,"Source_File"] %in% GenomesToAlign[i],1]))
       print(GenomesToAlign[i])
       print("")
       }
      }	
     }
    }
  }

  return(SamplesFQToGrab)
}

RunSSfqfetchPipeline <- function(SampleToGrab=BamsToRealign,SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipeLineLocations,Config="Config"){

  Pipeline <- getPipelinesPath("fqfetchpipeline")
  workflowExec <- getWorkflowParam("Executable")
  javaExec <- getExecPath("java")
  mode <- tolower(getWorkflowParam("Mode"))
  PipelineBase <- GetPipelinebase()
  rsyncExec <- getExecPath("rsync")
  FQDir <- PLs@FQDir
  bamDir <- PLs@BamDir  
  
  
  wgetExec <- "wget"
  CMDs_RSYNC <- "-av"
  CMDs_WGET <- ""
  CMDs2_RSYNC <- ""
  CMDs2_WGET <- "-P"    
  logFileCMD_WGET <- "-o "
  logFileCMD_RSYNC <- "--log-file="
  
  
  
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")
  

  SampleSheet[SampleSheet[,"BamLocation"] %in% c("/","NA"," "),"BamLocation"] <- "Location_Not_Known" 
  SampleSheet[SampleSheet[,"FQLocation"] %in% c("/","NA"," "),"FQLocation"] <- "Location_Not_Known" 
  
  genome <- GetGenomeFromConfig(WkgDir,ConfigDirectory=Config)
  Variables  <- c(file.path(WkgDir,""),FQDir,bamDir)
  names(Variables) <- c("WorkingDirectory","FQDir","bamDir")

  
  ToGrab <- SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe" & SampleSheet[,"GenomicsID"] %in% SampleToGrab & !SampleSheet[,"FQLocation"] %in% "Location_Not_Known","FQLocation"]
  NameToGrab <- SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe" & SampleSheet[,"GenomicsID"] %in% SampleToGrab  & !SampleSheet[,"FQLocation"] %in% "Location_Not_Known","GenomicsID"]
  
  FindThoseToBePlexed <- SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe" & SampleSheet[,"GenomicsID"] %in% unique(as.vector(SampleSheet[,"Source_File"])),"FQLocation"] 
  FindThoseNamesOfBeingPlexed <- SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe" & SampleSheet[,"GenomicsID"] %in% unique(as.vector(SampleSheet[,"Source_File"])),"GenomicsID"]
  
  ThoseWithFQNoBam <- SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe" & SampleSheet[,"BamLocation"] %in% "Location_Not_Known" & !SampleSheet[,"FQLocation"] %in% "Location_Not_Known","FQLocation"]
  NamesOfThoseWithFQNoBam <- SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe" & SampleSheet[,"BamLocation"] %in% "Location_Not_Known" & !SampleSheet[,"FQLocation"] %in% "Location_Not_Known","GenomicsID"]
  
  FindThoseNamesOfBeingPlexed <-   FindThoseNamesOfBeingPlexed[!is.na(FindThoseToBePlexed)]
  FindThoseToBePlexed <-   FindThoseToBePlexed[!is.na(FindThoseToBePlexed)]
  
  ThoseWithFQNoBam <-   ThoseWithFQNoBam[!is.na(ThoseWithFQNoBam)]
  NamesOfThoseWithFQNoBam <-   NamesOfThoseWithFQNoBam[!is.na(ThoseWithFQNoBam)]
  
  ToGrab <- c(ToGrab,FindThoseToBePlexed,ThoseWithFQNoBam)
  NameToGrab <- c(NameToGrab,FindThoseNamesOfBeingPlexed,NamesOfThoseWithFQNoBam)
  
  if(length(ToGrab) > 0){
    
    Specialisations  <- vector("list",length=length(ToGrab))
    for(i in 1:length(Specialisations)){
      if(grepl("^ftp://",as.vector(ToGrab[i]))){
        Specialisations[[i]] <- c(as.vector(ToGrab[i]),as.vector(NameToGrab[i]),wgetExec,CMDs_WGET,CMDs2_WGET,logFileCMD_WGET)
        names(Specialisations[[i]]) <- c("BamToGet","Bam_Name","rsync","CMDs","CMDs2","logFileCMD")        
      }
      if(!grepl("^ftp://",as.vector(ToGrab[i]))){
        Specialisations[[i]] <- c(as.vector(ToGrab[i]),as.vector(NameToGrab[i]),rsyncExec,CMDs_RSYNC,CMDs2_RSYNC,logFileCMD_RSYNC)
        names(Specialisations[[i]]) <- c("BamToGet","Bam_Name","rsync","CMDs","CMDs2","logFileCMD")        
      }
    } 
    names(Specialisations) <-  paste(as.vector(gsub("\\.bwa.*","",NameToGrab)),"_FQFetch",sep="") 
  
    #Pipeline = "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/FQGettingPipeline.xml"  
    PipeName <- "SS_FQFetch"
    SSFQFetch_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=MaxJobs)
    saveXML(SSFQFetch_WfMeta,file=file.path(PLs@WorkFlowDir,"SSFQFetch.xml"))
    system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,"SSFQFetch.xml"),sep=""),wait=T)
  }
  FQFilesDirs <- list.dirs(path=PLs@FQDir,recursive=F)
  BaseFQDirs  <-  basename(FQFilesDirs)
  for(i in 1:nrow(SampleSheet)){
  	if(!any(c(grep(".fastq$",SampleSheet[i,"Source_File"]),grep(".fastq.gz$",SampleSheet[i,"Source_File"]),grep(".fq.gz$",SampleSheet[i,"Source_File"]),grep(".fq$",SampleSheet[i,"Source_File"]),grep("Processed",SampleSheet[i,"Source_File"]),grep("Realign",SampleSheet[i,"Source_File"])))){
      if(any(NameToGrab %in% SampleSheet[i,"GenomicsID"])){
        if(any(BaseFQDirs %in% SampleSheet[i,"GenomicsID"])){
#  			SampleSheet[i,"Source_File"] <- SampleSheet[i,"FQLocation"][basename(SampleSheet[i,"FQLocation"]) %in% basename(dir(FQFilesDirs[BaseFQDirs %in% SampleSheet[i,"GenomicsID"]],full.names=T))]
        TempSampleFQbases <- basename(dir(FQFilesDirs[BaseFQDirs %in% SampleSheet[i,"GenomicsID"]],full.names=T))
        if(length(dir(FQFilesDirs[BaseFQDirs %in% SampleSheet[i,"GenomicsID"]],full.names=T)) > 0){
          SampleSheet[i,"Source_File"] <- dir(FQFilesDirs[BaseFQDirs %in% SampleSheet[i,"GenomicsID"]],full.names=T)[TempSampleFQbases %in% basename(SampleSheet[i,"FQLocation"])]
        }
  		}else{
  		    SampleSheet[i,"Source_File"] <- "NA"
  		    SampleSheet[i,"Analysis_State"] <- "FQ required!!"
  		}	
  	 }
    }
  }
  return(SampleSheet)  
}
RunSSdeMultiplesPipeline <- function(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipeLineLocations,Config="Config"){

    DemultiplexSS <- SampleSheet[SampleSheet[,"BamLocation"] %in% "Premultiplex",]
    AllSourceFiles <- as.vector(DemultiplexSS[,"Source_File"])    
    dir.create(file.path(PLs@FQDir),showWarnings = FALSE)
    Specialisations  <- vector("list")
    if(nrow(DemultiplexSS) > 0){

    for(i in 1:length(unique(DemultiplexSS[,"Source_File"]))){
      FQwhichBarcodesFrom <- as.vector(SampleSheet[SampleSheet[,"GenomicsID"] %in% unique(DemultiplexSS[,"Source_File"])[i],"Source_File"])
      SequenceExpected <-  as.vector(SampleSheet[SampleSheet[,"GenomicsID"] %in% unique(DemultiplexSS[,"Source_File"])[i],"GenomicsID"])
      FQDir <- file.path(PLs@FQDir,as.vector(SampleSheet[SampleSheet[,"GenomicsID"] %in% unique(DemultiplexSS[,"Source_File"])[i],"GenomicsID"]),"")
      print(FQwhichBarcodesFrom)
      print(SequenceExpected)
      print(FQDir)
      if(length(FQwhichBarcodesFrom) > 0 & length(SequenceExpected) > 0 & length(FQDir) > 0){    
        if(grep("^SLX",SequenceExpected) > 0){
          BarcodeFile <- dir(path=file.path(PLs@LocationsDir,"barcodes"),full.names=T)[gsub(".txt$","",dir(path=file.path(PLs@LocationsDir,"barcodes"),full.names=F)) %in% gsub("\\..*","",SequenceExpected)]
        }else{
          BarcodeFile <- dir(path=file.path(PLs@LocationsDir,"barcodes"),full.names=T)[gsub(".txt$","",dir(path=file.path(PLs@LocationsDir,"barcodes"),full.names=F)) %in% SequenceExpected]
     
        }      
  
        Specialisations <- c(Specialisations,list(c(FQwhichBarcodesFrom,SequenceExpected,FQDir,BarcodeFile)))
        names(Specialisations[[length(Specialisations)]]) <- c("SequenceFile","SequenceExpected","demuxDir","indexFile")
      }
      
      }
      names(Specialisations) <- seq(1,length(Specialisations))
      Variables <- WkgDir
      names(Variables) <- "WorkingDirectory"
      Pipeline <- "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/SpecialisationDemultiplex.xml"
      PipeName <- "Demultiplex"
      SSDemultiplex_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=75)
      saveXML(SSDemultiplex_WfMeta,file=file.path(PLs@WorkFlowDir,"SSDemultiplex.xml"))
      system(paste("java -jar /lustre/mib-cri/carrol09/MyPipe/workflow-all-1.2-SNAPSHOT.jar --mode=lsf ",file.path(PLs@WorkFlowDir,"SSDemultiplex.xml"),sep=""),wait=T) 
      }       
      for(i in 1:nrow(SampleSheet)){
        if(SampleSheet[i,"FQLocation"] %in% "Premultiplex"){
          AllFileForSample <- dir(path=file.path(PLs@FQDir,SampleSheet[i,"Source_File"]),pattern=".gz$")
          AllFileForSampleFull <- dir(path=file.path(PLs@FQDir,SampleSheet[i,"Source_File"]),pattern=".gz$",full.names=T)
          SampleSheet[i,"FQLocation"] <- AllFileForSampleFull[gsub(".gz$","",AllFileForSample) %in% SampleSheet[i,"SampleName"]]
          SampleSheet[i,"Source_File"] <- AllFileForSampleFull[gsub(".gz$","",AllFileForSample) %in% SampleSheet[i,"SampleName"]]          
          SampleSheet[i,"BamLocation"] <- "Pre-Alignment"
        }
      }
      return(SampleSheet)
}

RunSSRealignmentPipeline <- function(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipeLineLocations,Config="Config"){


  Pipeline <- getPipelinesPath("alignpipeline")
  workflowExec <- getWorkflowParam("Executable")
  javaExec <- getExecPath("java")
  mode <- tolower(getWorkflowParam("Mode"))
  PipelineBase <- GetPipelinebase()
  bwaExec <- getExecPath("bwa")
  pythonExec <- getExecPath("python")
  samtoolsExec <- getExecPath("samtools")  
  picardExec <- getExecPath("picard")    
  pythonlibs <- getLibraryPath("pythonlibs")
  fastxtrimmerexe <- getFastxTrimmerExec()
  maxlength <- getReadLengthMax()
  #QE <- GetQualityEncodingFromConfig()
  QEflag <- "NoValueNeeded"


  FQDir <- PLs@FQDir
  bamDir <- PLs@BamDir  

  
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")

  RemainderSheet <- SampleSheet[!SampleSheet[,"Analysis_State"] %in% "RunMe",]
  SampleSheet <- SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe",]

  if(length(grep("fq$|.gz$|.gz$",SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe","Source_File"])) > 1){   
    SampleSheetOfThoseToAlign <- SampleSheet[grep("fq$|.gz$|.gz$",SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe","Source_File"]),]
    SampleSheetOfThoseToAlign <- SampleSheetOfThoseToAlign[!is.na(SampleSheetOfThoseToAlign[,"Source_File"]) & !SampleSheetOfThoseToAlign[,"Source_File"] %in% "NA",]
  }else{
    SampleSheetOfThoseToAlign <- matrix(SampleSheet[grep("fq$|.gz$|.gz$",SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe","Source_File"]),],ncol=ncol(SampleSheet))
    colnames(SampleSheetOfThoseToAlign) <- colnames(SampleSheet)
    SampleSheetOfThoseToAlign <- matrix(SampleSheetOfThoseToAlign[!is.na(SampleSheetOfThoseToAlign[,"Source_File"]) & !SampleSheetOfThoseToAlign[,"Source_File"] %in% "NA",],ncol=ncol(SampleSheet))
    colnames(SampleSheetOfThoseToAlign) <- colnames(SampleSheet)    
  }



  TargetGenome <- GetGenomeFromConfig(WkgDir,ConfigDirectory="Config")
  GenomeBuild <- GetGenomeBuildFromConfig(WkgDir,ConfigDirectory="Config")
  Specialisations  <- vector("list")
  
  Variables <- c(WkgDir,FQDir,TargetGenome,GenomeBuild,bwaExec,pythonExec,javaExec,samtoolsExec,PipelineBase,picardExec,pythonlibs,bamDir,fastxtrimmerexe,maxlength,QEflag)
  names(Variables) <- c("WorkingDirectory","FQDirectory","GenomeBuild","Genome","bwa","python","java","samtools","pipelineBase","picard","pythonlibs","bamDir","fastxtrimmerexe","maxlength","QEflag")
  if(nrow(SampleSheetOfThoseToAlign) > 0){ 
    for(i in 1:nrow(SampleSheetOfThoseToAlign)){
      SampleID <- SampleSheetOfThoseToAlign[i,"GenomicsID"]
      FileToAlign <- SampleSheetOfThoseToAlign[i,"Source_File"]
      Specialisations <- c(Specialisations,list(c(SampleID,FileToAlign)))
      names(Specialisations[[length(Specialisations)]]) <- c("Test","FastaFile")  
    }
     names(Specialisations) <- seq(1,length(Specialisations))  
#     Pipeline <- "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/ReAlignPipeline.xml"
     PipeName <- "Alignment"
     SSRealignment_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=MaxJobs)
     saveXML(SSRealignment_WfMeta,file=file.path(PLs@WorkFlowDir,"SSAlignment.xml"))
     system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,"SSAlignment.xml"),sep=""),wait=T)
   }
   if(nrow(SampleSheet) > 0){
     for(i in 1:nrow(SampleSheet)){
       if(file.exists(file.path(PLs@BamDir,paste(SampleSheet[i,"GenomicsID"],".bwa.Realigned",TargetGenome,".bam",sep="")))){
          SampleSheet[i,"Source_File"] <- paste(SampleSheet[i,"GenomicsID"],".bwa.Realigned",TargetGenome,".bam",sep="") 
       }
     }
        SampleSheet <- rbind(SampleSheet,RemainderSheet)
    }else{
        SampleSheet <- RemainderSheet
    }
   return(SampleSheet)
     
}


RunSSMergingPipeline <- function(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipeLineLocations,Config="Config"){

  Pipeline <- getPipelinesPath("mergingpipeline")
  workflowExec <- getWorkflowParam("Executable")
  javaExec <- getExecPath("java")
  mode <- tolower(getWorkflowParam("Mode"))
  PipelineBase <- GetPipelinebase()
  picardExec <- getExecPath("picard")    
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")
  bamDir <- PLs@BamDir

    TempSampleSheet <- SampleSheet
    UniquedMergedFiles <-  as.vector(unique(na.omit(SampleSheet[,"toMerge"])))
    MoreThan1Sample <- names(table(na.omit(SampleSheet[,"SampleName"])))[table(na.omit(SampleSheet[,"SampleName"])) > 1]
    AllSamplesPresent <- as.vector(unique(na.omit(SampleSheet[,"GenomicsID"])))


 
    
    UniquedMergedFiles <- UniquedMergedFiles[!UniquedMergedFiles %in% c(AllSamplesPresent)]
    UniquedMergedFiles <- UniquedMergedFiles[!UniquedMergedFiles %in% MoreThan1Sample]
         
    SampleSheet <- SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe",,drop=F]     
    InputLists <- vector("character")
    InputLists2 <- vector("character")
    OutputNames <- vector("character")          
    OutputNames2 <- vector("character")    
    if(length(UniquedMergedFiles) != 0){
 
      OutputNames <- UniquedMergedFiles
      for(i in 1:length(UniquedMergedFiles)){
          if(length(grep("*.bam$",as.vector(SampleSheet[SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe","toMerge"] %in% UniquedMergedFiles[i],"Source_File"]))) > 1){
            ToBeMerged <-  as.vector(SampleSheet[SampleSheet[,"toMerge"] %in% UniquedMergedFiles[i],"Source_File"])
            ToBeMerged <- ToBeMerged[!is.na(ToBeMerged)]
            ToBeMerged <- ToBeMerged[!grepl(".fq$|fq.gz$",ToBeMerged)]
            if(length(ToBeMerged) > 0){
              InputLists[i] <- paste(paste("INPUT=",file.path(PLs@BamDir,ToBeMerged),sep=""),collapse=" ")
            }
          }
      }
      UniquedMergedFiles <- UniquedMergedFiles[!gsub("INPUT=","",InputLists) == ""]   
    }
    
    MoreThan1Sample <- MoreThan1Sample[!MoreThan1Sample %in% c(AllSamplesPresent)]  
    if(AutoMergeCheck(WkgDir,"Config")){
    if(length(MoreThan1Sample) != 0){
      InputLists2 <- vector("character",length=length(MoreThan1Sample))
      OutputNames2 <- MoreThan1Sample
      for(i in 1:length(MoreThan1Sample)){
        if(length(grep("*.bam$",as.vector(SampleSheet[SampleSheet[,"SampleName"] %in% MoreThan1Sample[i],"Source_File"]))) > 1){
          ToBeMerged2 <-  as.vector(SampleSheet[SampleSheet[,"SampleName"] %in% MoreThan1Sample[i],"Source_File"])
          ToBeMerged2 <- ToBeMerged2[!is.na(ToBeMerged2)]
          ToBeMerged2 <- ToBeMerged2[!grepl(".fq$|fq.gz$",ToBeMerged2)]
                    ToBeMerged2 <- ToBeMerged2[!grepl("Not_Found",ToBeMerged2)] 
          if(length(ToBeMerged2) > 0){
             InputLists2[i] <- paste(paste("INPUT=",file.path(PLs@BamDir,ToBeMerged2),sep=""),collapse=" ")
          }         
        }
      }
      MoreThan1Sample <- MoreThan1Sample[!gsub("INPUT=","",InputLists2) == ""]
    }
    InputLists <- c(InputLists,InputLists2)
    OutputNames <- c(OutputNames,OutputNames2)
    }
    Variables <- c(WkgDir,picardExec,javaExec,bamDir)
    names(Variables) <- c("WorkingDirectory","picard","java","bamDir")
    Specialisations  <- vector("list")
  if(length(InputLists) > 0){       
    for(i in 1:length(InputLists)){   
        Inputs <- InputLists[i]
         if(Inputs != ""){
        Outputs <- OutputNames[i]
        Specialisations <- c(Specialisations,list(c(Inputs,Outputs)))
        names(Specialisations[[length(Specialisations)]]) <- c("InputName","OutName")  
      }
    }
    if(length(Specialisations) > 0){
     names(Specialisations) <- seq(1,length(Specialisations))  
#     Pipeline <- "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/MergingPipeline.xml"
     PipeName <- "Merging"
     SSMerge_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=75)
     saveXML(SSMerge_WfMeta,file=file.path(PLs@WorkFlowDir,"SSMerging.xml"))
     system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,"SSMerging.xml"),sep=""),wait=T)
     SampleSheet <- TempSampleSheet
     for(i in 1:length(OutputNames)){
     TempMatrix <- matrix(nrow=length(OutputNames),ncol=ncol(SampleSheet),data=NA)
     colnames(TempMatrix) <- colnames(SampleSheet)
         if(file.exists(file.path(PLs@BamDir,paste(OutputNames[i],".bwa.bam",sep="")))){
            TempMatrix[i,"GenomicsID"] <- OutputNames[i]
            TempMatrix[i,"SampleName"] <- OutputNames[i]
            TempMatrix[i,"BamLocation"] <- file.path(PLs@BamDir,paste(OutputNames[i],".bwa.bam",sep=""))
            TempMatrix[i,"Source_File"] <- paste(OutputNames[i],".bwa.bam",sep="") 
            TempMatrix[i,"Analysis_State"]  <- "RunMe"            
            SampleSheet <- rbind(SampleSheet,TempMatrix)
         }
     }
     SampleSheet[SampleSheet[,"toMerge"] %in% OutputNames ,"Analysis_State"]  <- "MergedAndUnused"
     SampleSheet[!SampleSheet[,"GenomicsID"] %in% OutputNames & SampleSheet[,"SampleName"] %in% OutputNames & SampleSheet[,"Analysis_State"] %in% "RunMe","Analysis_State"]  <- "MergedAndUnused"
   }
   }
   SampleSheet <- SampleSheet[!is.na(SampleSheet[,"GenomicsID"]),]
   SampleSheet
}

RunBamProcessPipeline <- function(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipelineLocations,Config="Config"){
   
#  genome <- GetGenomeFromConfig(WkgDir,ConfigDirectory=Config)

  Pipeline <- getPipelinesPath("bamprocesspipeline")
  workflowExec <- getWorkflowParam("Executable")
  javaExec <- getExecPath("java")
  mode <- tolower(getWorkflowParam("Mode"))
  PipelineBase <- GetPipelinebase()
  pythonExec <- getExecPath("python")    
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")
  rExec <- getExecPath("rexec")
  pythonlibs <- getLibraryPath("pythonlibs")
  rlibs <- getLibraryPath("rlibs")
  FLDir <- PLs@FragLengthDir



  MapQFlag <- GetMapQFlag(WkgDir,Config)
  DupFlag <- GetDupFlag(WkgDir,Config)
  ExcludedFlag <- GetExcludedFlag(WkgDir,Config)  
  
  Variables  <- c(file.path(WkgDir,""),PLs@BamDir,"Dummy",DupFlag,ExcludedFlag,MapQFlag,pythonExec,PipelineBase,rExec,pythonlibs,rlibs,FLDir)
  names(Variables) <- c("WorkingDirectory","BamDirectory","genomeFile","DupFlag","ExcludedFlag","MapQFlag","python","pipelineBase","rexec","pythonlibs","rlibs","FLDir")
  
  BamFiles <- SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe" & !SampleSheet[,"Source_File"] %in% "NA" & !is.na(SampleSheet[,"Source_File"]) & grepl(".bam",SampleSheet[,"Source_File"]) & !grepl("_Processed.bam",SampleSheet[,"Source_File"]),"Source_File"]
  BamFiles <- gsub(".bam","",BamFiles)
  
  genome <- GetGenomeFromConfig(WkgDir,ConfigDirectory="Config")
  excludedList <- GetExcludedListFromConfig(WkgDir,ConfigDirectory="Config")
  SeqDict <- GetSequenceDictionaryFromConfig(WkgDir,ConfigDirectory="Config")
  
  if(length(BamFiles) > 0){
    Specialisations  <- vector("list",length=length(BamFiles))
    for(i in 1:length(Specialisations)){
      Specialisations[[i]] <- c(as.vector(BamFiles[i]),genome,excludedList,SeqDict)
      names(Specialisations[[i]]) <- c("Test","genome","excludedfile","seqdict")  
    } 
    names(Specialisations) <-  paste(as.vector(gsub("\\.bwa.*","",BamFiles)),"_bamProcess",sep="") 
    #Specialisations[[1]] <- "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/20121114_MontoyaAR_DN_Hes6ChIP/bamFiles/ID-LNCAP-LM-HES6-BICALUTAMIDE-INPUT-D26.bwa.bam"
  
#    Pipeline = "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/BamProcessPipeline.xml"  
    PipeName <- "SS_BamProcess"
    SSBamProcess_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=MaxJobs)
    saveXML(SSBamProcess_WfMeta,file=file.path(PLs@WorkFlowDir,"SS_BamProcess.xml"))
    system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,"SS_BamProcess.xml"),sep=""),wait=T)
  }
  fileLogs <- dir(path=PLs@BamDir,pattern="_fileLog.log",full.names = F)
  fileLogsFull <- dir(path=PLs@BamDir,pattern="_fileLog.log",full.names = T)  
  
  CovAndSisserLogs <- dir(path=PLs@FragLengthDir,pattern="AllFragLog$",full.names = F)
  CovAndSisserLogsFull <- dir(path=PLs@FragLengthDir,pattern="AllFragLog$",full.names = T)

  CorrLogs <- dir(path=PLs@FragLengthDir,pattern="CorrFragLog$",full.names = F)
  CorrLogsFull <- dir(path=PLs@FragLengthDir,pattern="CorrFragLog$",full.names = T)  
  
  
  BamsInSampleSheet <-  as.vector(SampleSheet[,"Source_File"])
  SampleSheet <- data.frame(SampleSheet,stringsAsFactors=F)
  
  for(i in 1:nrow(SampleSheet)){
    if(any(gsub("_fileLog.log","",fileLogs) %in% gsub(".bam","",SampleSheet[i,"Source_File"]))){
       LogToRead <- fileLogsFull[gsub("_fileLog.log","",fileLogs) %in% gsub(".bam","",SampleSheet[i,"Source_File"])]
       DataIn <- read.delim(LogToRead,sep="\t")
       SampleSheet[i,c("Original","delRand","Excluded","Filtered","Unique")] <- DataIn[,c("Mapped","NonRandomChr","IncludedRegions","QC...15","Unique")]
       SampleSheet[i,"DuplicationRate"] <- (as.numeric(SampleSheet[i,"Filtered"])-as.numeric(SampleSheet[i,"Unique"]))/as.numeric(SampleSheet[i,"Filtered"])*100
       SampleSheet[i,"NRF"] <- DataIn[,"Unique"]/DataIn[,"QC...15"]

    }else{
    #  SampleSheet[i,"Source_File"] <- "No_Processed_Bam"
    	SampleSheet[i,c("Original","delRand","Excluded","Filtered","Unique")] <- rep("No_Information_Available",5)
    	SampleSheet[i,"DuplicationRate"] <- "No_Information_Available"
    }
    if(any(gsub("_Processed.AllFragLog","",CovAndSisserLogs) %in% gsub(".bam","",SampleSheet[i,"Source_File"]))){
    		FragsToRead <- CovAndSisserLogsFull[gsub("_Processed.AllFragLog","",CovAndSisserLogs) %in% gsub(".bam","",SampleSheet[i,"Source_File"])]
    		print(FragsToRead)
    		DataIn <- as.matrix(read.delim(FragsToRead,sep=" ",comment.char="#"))
    		SampleSheet[i,"Sissr_Fragment_Length"] <- DataIn[1,1]		
    		SampleSheet[i,"Correlation_Fragment_Length"] <- DataIn[1,2]
    		SampleSheet[i,"Coverage_Fragment_Length"] <- DataIn[1,3]
   	}else{
    		SampleSheet[i,"Sissr_Fragment_Length"] <- "No_Information_Available"		
    		SampleSheet[i,"Correlation_Fragment_Length"] <- "No_Information_Available"
    		SampleSheet[i,"Coverage_Fragment_Length"] <- "No_Information_Available"	
    
   	}
    if(any(gsub("_Processed.CorrFragLog","",CorrLogs) %in% gsub(".bam","",SampleSheet[i,"Source_File"]))){
    		CorrToRead <- CorrLogsFull[gsub("_Processed.CorrFragLog","",CorrLogs) %in% gsub(".bam","",SampleSheet[i,"Source_File"])]
    		print(CorrToRead)
    		DataIn <- as.matrix(read.delim(CorrToRead,sep=" ",comment.char="#"))
    		SampleSheet[i,"CC_Fragment_Length"] <- DataIn[1,1]		
    		SampleSheet[i,"NSC"] <- DataIn[1,2]
    		SampleSheet[i,"RSC"] <- DataIn[1,3]
   	}else{
    		SampleSheet[i,"CC_Fragment_Length"] <- "No_Information_Available"		
    		SampleSheet[i,"NSC"] <- "No_Information_Available"
    		SampleSheet[i,"RSC"] <- "No_Information_Available"	
   	}
   	if(any(gsub("_fileLog.log","",fileLogs) %in% gsub(".bam","",SampleSheet[i,"Source_File"]))){
      if(file.exists(gsub("_fileLog.log$","_Processed.bam",fileLogsFull)[gsub("_fileLog.log$","",fileLogs) %in% gsub(".bam","",SampleSheet[i,"Source_File"])])){
            SampleSheet[i,"Processed_bamFileName"] <- gsub("_fileLog.log$","_Processed.bam",fileLogs)[gsub("_fileLog.log$","",fileLogs) %in% gsub(".bam","",SampleSheet[i,"Source_File"])]
      }else{
                  # SampleSheet[i,"Source_File"] <- "No_Processed_Bam"   
      }
    }
  }
  return(SampleSheet)
}



RunBamProfilePipeline <- function(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipelineLocations,Config="Config"){


  Pipeline <- getPipelinesPath("bamprofilepipeline")
  workflowExec <- getWorkflowParam("Executable")
  javaExec <- getExecPath("java")
  mode <- tolower(getWorkflowParam("Mode"))
  PipelineBase <- GetPipelinebase()
  pythonExec <- getExecPath("python")    
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")
  rexec <- getExecPath("rexec")
  bedtools <- getExecPath("bedtools")
  BigWig <- getExecPath("bigwig")
  pythonlibs <- getLibraryPath("pythonlibs")
  rlibs <- getLibraryPath("rlibs")


  CoverageDir <- PLs@CovDir
  GenePos <- GetGenePosFromConfig(WkgDir,ConfigDirectory="Config")
  GeneSets <- GetGeneSetsFromConfig(WkgDir,ConfigDirectory="Config")


   
   genome <- GetGenomeFromConfig(WkgDir,ConfigDirectory="Config") 
   genomeChrLengths <- GetChrLengthsFromConfig(WkgDir,ConfigDirectory="Config")
   
#  genome <- GetGenomeFromConfig(WkgDir,ConfigDirectory=Config)
  Variables  <- c(file.path(WkgDir,""),PLs@BamDir,genome,genomeChrLengths,BigWig,bedtools,rexec,PipelineBase,pythonlibs,rlibs,GenePos,GeneSets,CoverageDir)
  names(Variables) <- c("WorkingDirectory","BamDirectory","genomeName","genomeFile","BigWig","bedtools","rexec","pipelineBase","pythonlibs","rlibs","genepos","genesets","CoverageDir")
  
  BamProfileFiles <- SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe" & SampleSheet[,"Processed_bamFileName"] != "NA" & !is.na(SampleSheet[,"Processed_bamFileName"]) & grepl(".bam",SampleSheet[,"Processed_bamFileName"]),"Processed_bamFileName"]
  BamProfileFiles <- gsub(".bam","",BamProfileFiles)
  
  Specialisations  <- vector("list",length=length(BamProfileFiles))
  for(i in 1:length(Specialisations)){
    Specialisations[[i]] <- c(as.vector(BamProfileFiles[i]))
    names(Specialisations[[i]]) <- c("Test")  
  } 
  names(Specialisations) <-  paste(as.vector(gsub("\\.bwa.*","",BamProfileFiles)),"_bamProcess",sep="") 

#  Pipeline = "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/MultiBamProcessPipeline_P2.xml"  
  PipeName <- "SS_BamProfile"
  SSBamProfile_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=MaxJobs)
  saveXML(SSBamProfile_WfMeta,file=file.path(PLs@WorkFlowDir,"SS_BamProfile.xml"))
  system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,"SS_BamProfile.xml"),sep=""),wait=T)


    SampleSheet <- ReadAndLock(file.path(WkgDir,"SampleSheet.csv"),WkdDir,SAF=T,napTime=5)
    
    
  tryCatch({ 
    SampleSheetTemp <- SampleSheet
    
    for(i in 1:nrow(SampleSheetTemp)){
    	SampleToLookFor <- gsub(".bam","",SampleSheetTemp[i,"Processed_bamFileName"])
    	BedGraphFile <- dir(path=CoverageDir,pattern=paste(SampleToLookFor,".bedgraph",sep=""),full.names=T)
    	if(length(BedGraphFile) > 0){
    		SampleSheetTemp[i,"BedGraph_Files"] <- BedGraphFile 
    	}
    
    	BigWigFile <- dir(path=CoverageDir,pattern=paste(SampleToLookFor,".bw",sep=""),full.names=T)
    	if(length(BigWigFile) > 0){
    		SampleSheetTemp[i,"BigWig_Files"] <- BigWigFile
    	}
    
    	CovStats <- dir(path=CoverageDir,pattern=paste(SampleToLookFor,"_CoverageDispersion.RData",sep=""),full.names=T)
    	if(length(CovStats) > 0){
    		load(CovStats)
    		SampleSheetTemp[i,"SSD_Of_Coverage"] <- ssdOfSample
    		SampleSheetTemp[i,"Gini_Of_Coverage"] <- giniOfSample[1]
    		SampleSheetTemp[i,"adjusted_Gini_Of_Coverage"] <- giniOfSample[2]
    		rm(giniOfSample)
    		rm(ssdOfSample)
    	}
    	ReadsInFeatures <- dir(path=CoverageDir,pattern=paste(SampleToLookFor,"_ReadCountsInFeatures.txt",sep=""),full.names=T)
    	if(length(ReadsInFeatures) > 0){
    	
    		Temp <- read.delim(ReadsInFeatures,h=T,sep="\t")
    		if(any(colnames(Temp) %in% "TSS_500Upstream_500Downstream")){
    			Total <- Temp["Counts","Total"]
    			TSS <- Temp["Counts","TSS_500Upstream_500Downstream"]
    			Pro2000_500 <- Temp["Counts","Promoter_2000Upstream_500Upstream"]
    			Pro5000_2000 <- Temp["Counts","Promoter_5000Upstream_2000Upstream"]
    			Pro10000_5000 <- Temp["Counts","Promoter_10000Upstream_5000Upstream"]
    			GeneBody <- Temp["Counts","GeneBody_minus_TSS"]
    			Intergenic <- Total-(TSS + Pro2000_500 + Pro5000_2000 + Pro10000_5000 + GeneBody)
    			TSSPer <- (TSS/Total)*100
    			Pro2000_500Per <- (Pro2000_500/Total)*100
    			Pro5000_2000Per <- (Pro5000_2000/Total)*100
    			Pro10000_5000Per <- (Pro10000_5000/Total)*100
    			GeneBodyPer <- (GeneBody/Total)*100
    			IntergenicPer <- (Intergenic/Total)*100
    SampleSheetTemp[i,c("Reads_In_TSS")] <- TSS
    SampleSheetTemp[i,c("Reads_In_Promoter_2000Up_500Up")] <- Pro2000_500
    SampleSheetTemp[i,c("Reads_In_Promoter_5000Up_2000Up")] <- Pro5000_2000
    SampleSheetTemp[i,c("Reads_In_Promoter_10000Up_5000Up")] <- Pro10000_5000
    SampleSheetTemp[i,c("Reads_In_GeneBody_Minus_TSS")] <- GeneBody
    SampleSheetTemp[i,c("Reads_In_Intergenic")]	<- Intergenic	
    
    SampleSheetTemp[i,c("Percent_Of_Reads_In_TSS")] <- TSSPer
    SampleSheetTemp[i,c("Percent_Of_Reads_In_Promoter_2000Up_500Up")] <- Pro2000_500Per
    SampleSheetTemp[i,c("Percent_Of_Reads_In_Promoter_5000Up_2000Up")] <- Pro5000_2000Per
    SampleSheetTemp[i,c("Percent_Of_Reads_In_Promoter_10000Up_5000Up")] <- Pro10000_5000Per
    SampleSheetTemp[i,c("Percent_Of_Reads_In_GeneBody_Minus_TSS")] <- GeneBodyPer
    SampleSheetTemp[i,c("Percent_Of_Reads_In_Intergenic")]	<- IntergenicPer			
    		}
    	}
    
    }	
    
   SampleSheet <- SampleSheetTemp 	
  },warning = function(war){
    write.table(war,file.path(PLs@WorkFlowDir,paste(PipeName,"_WarningFile.txt",sep="")),sep="\t")
  
  },error = function(err){
    write.table(err,file.path(PLs@WorkFlowDir,paste(PipeName,"_ErrorFile.txt",sep="")),sep="\t")
  
  },finally = {
    WriteAndUnlock(SampleSheet,file.path(WkgDir,"SampleSheet.csv"))
  }
  
  )
  
}

RunMacsPeakCallerPipeline <- function(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipelineLocations,Config="Config"){
if(CallPeaksCheck("Macs",WkgDir,Config)){
   
   Mfold <- getMacsmfold(WkgDir,ConfigDirectory="Config")
   MacsGenome = getMacsGenome(WkgDir,ConfigDirectory="Config")
   ShiftSizeDefault = getMacsShiftSize(WkgDir,ConfigDirectory="Config")    

  Pipeline <- getPipelinesPath("macspeakcallpipeline")
  workflowExec <- getWorkflowParam("Executable")
  javaExec <- getExecPath("java")
  mode <- tolower(getWorkflowParam("Mode"))
  PipelineBase <- GetPipelinebase()
  pythonExec <- getExecPath("python")    
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")
  rexec <- getExecPath("rexec")
  bedtools <- getExecPath("bedtools")
  BigWig <- getExecPath("bigwig")
  macs <- getExecPath("macs")
  pythonlibs <- getLibraryPath("pythonlibs")

   

   
   genome <- GetGenomeFromConfig(WkgDir,ConfigDirectory="Config") 
   #genomeChrLengths <- GetChrLengthsFromConfig(WkgDir,ConfigDirectory="Config")
   SampleSheet <- SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe",]
   SamplesAndInputs <- matrix(data=c(gsub("\\.bam","",SampleSheet[,"Processed_bamFileName"]),gsub("\\.bam","",SampleSheet[match(SampleSheet[,"InputToUse"],SampleSheet[,"SampleName"],incomparable=NA),"Processed_bamFileName"]),as.vector(SampleSheet[,"CC_Fragment_Length"])),ncol=3,byrow=F,dimnames=list(NULL,c("Samples","Inputs","CC_Fraglength")))
   SamplesAndInputs <- SamplesAndInputs[!is.na(SamplesAndInputs[,"Inputs"]) & ! as.vector(SamplesAndInputs[,"Samples"]) %in% "No_Processed_Bam",]

   SamplesAndInputs2 <- matrix(data=c(gsub("\\.bam","",SampleSheet[,"Processed_bamFileName"]),gsub("\\.bam","",SampleSheet[match(SampleSheet[,"InputToUse"],SampleSheet[,"GenomicsID"],incomparable=NA),"Processed_bamFileName"]),as.vector(SampleSheet[,"CC_Fragment_Length"])),ncol=3,byrow=F,dimnames=list(NULL,c("Samples","Inputs","CC_Fraglength")))  
   SamplesAndInputs2 <- SamplesAndInputs2[!is.na(SamplesAndInputs2[,"Inputs"]) & ! as.vector(SamplesAndInputs2[,"Samples"]) %in% "No_Processed_Bam",]

   MakeUniquePeakCalls <- rbind(SamplesAndInputs,SamplesAndInputs2)
   MakeUniquePeakCalls <- MakeUniquePeakCalls[match(unique(MakeUniquePeakCalls[,1]),MakeUniquePeakCalls[,1]),]   
   SamplesAndInputs <- matrix(MakeUniquePeakCalls,ncol=3,byrow=F)

   SamplesAndInputs[SamplesAndInputs[,3] %in% "Too_few_Reads_To_Calculate" | SamplesAndInputs[,3] %in% "No_Information_Available" |  is.na(SamplesAndInputs[,3]),3] <-  ShiftSizeDefault

   for(i in 1:nrow(SamplesAndInputs)){
         if(!is.na(SampleSheet[gsub(".bam","",SampleSheet[,"Processed_bamFileName"]) %in% SamplesAndInputs[i,1],"InputFileNameToUse"])){
              SamplesAndInputs[i,2] <- gsub(".bam","",as.vector(SampleSheet[gsub(".bam","",SampleSheet[,"Processed_bamFileName"]) %in% SamplesAndInputs[i,1],"InputFileNameToUse"]))
         }
   }

  Variables  <- c(file.path(PLs@MacsDir,""),file.path(PLs@BamDir,""),MacsGenome,Mfold,macs,rexec,pythonlibs,pythonExec)
  names(Variables) <- c("Macs_Directory","BamDirectory","Genome","Mfold","macs","rexec","pythonlibs","python")
  

  if(nrow(SamplesAndInputs) > 0){
    Specialisations  <- vector("list",length=nrow(SamplesAndInputs))
    for(i in 1:length(Specialisations)){
      Specialisations[[i]] <- c(as.vector(SamplesAndInputs[i,1]),as.vector(SamplesAndInputs[i,2]),as.vector(SamplesAndInputs[i,3]))
      names(Specialisations[[i]]) <- c("Test","Control","shiftSize")  
    } 
    names(Specialisations) <-  paste(as.vector(SamplesAndInputs[,1]),"_MacsPeakCall",sep="") 
  
    #Pipeline = "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/MacsPeakCallingPipeline.xml"  
    PipeName <- "SS_MacsPeak"
    SSMacsPeak_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=75)
    saveXML(SSMacsPeak_WfMeta,file=file.path(PLs@WorkFlowDir,"SS_MacsPeak.xml"))
    system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,"SS_MacsPeak.xml"),sep=""),wait=T)
}


  MacsFiles <- dir(path=PLs@MacsDir,pattern="*peaks.xls$",full.name=F)
  MacsFilesFull <- dir(path=PLs@MacsDir,pattern="*peaks.xls$",full.name=T)  
  PositivePeaks <- MacsFiles[-grep("negative_peaks.xls",MacsFiles)]
  PositivePeaksFull <- MacsFilesFull[-grep("negative_peaks.xls",MacsFilesFull)]

  SampleSheet <- ReadAndLock(file.path(WkgDir,"SampleSheet.csv"),WkdDir,SAF=F,napTime=5)

      
  tryCatch({ 
    SampleSheetTemp <- SampleSheet
    

  for(i in 1:nrow(SampleSheetTemp)){
	 if(gsub(".bam","",SampleSheetTemp[i,"Processed_bamFileName"]) %in% gsub("_peaks.xls","",PositivePeaks)){
	   	MacsToRead <- file.path(PLs@MacsDir,PositivePeaks[gsub("_peaks.xls","",PositivePeaks) %in% gsub(".bam","",SampleSheetTemp[i,"Processed_bamFileName"])])
 	    print(MacsToRead)
		  DataIn <- read.delim(MacsToRead,sep="\t",comment.char="#")
		  SampleSheetTemp[i,"MacsPeaks"] <- nrow(DataIn)		
		  SampleSheetTemp[i,"Macs_name"] <- MacsToRead
		  print(nrow(DataIn))
		  if(file.exists(gsub("_peaks.xls","_model.r",MacsToRead))){
		    TempForFragLengths <- readLines(gsub("_peaks.xls","_model.r",MacsToRead))
		    SampleSheetTemp[i,"Macs_Fragment_Length"] <- gsub("d=","",gsub("\".*","",gsub(".*right\",c\\(\"d=","",TempForFragLengths[grep("d=",TempForFragLengths)])))
		  }else{
		     SampleSheetTemp[i,"Macs_Fragment_Length"] <- "Default:No_Model"
		     write.table("d=Default:No_Model",gsub("_peaks.xls","_model.r",MacsToRead),sep="",row.names=F,col.names=F,quote=T)
		     write.table("Dummy",gsub("_peaks.xls","_model.pdf",MacsToRead),sep="",row.names=F)		     
		  }
  	}
  }

   SampleSheet <- SampleSheetTemp 	
  },warning = function(war){
    write.table(war,file.path(PLs@WorkFlowDir,paste(PipeName,"_WarningFile.txt",sep="")),sep="\t")
  
  },error = function(err){
    write.table(err,file.path(PLs@WorkFlowDir,paste(PipeName,"_ErrorFile.txt",sep="")),sep="\t")
  
  },finally = {
    WriteAndUnlock(SampleSheet,file.path(WkgDir,"SampleSheet.csv"))
  }
  
  )

        return(SampleSheet)
}        
}

RunSicerPeakCallerPipeline <- function(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipelineLocations,Config="Config"){
if(CallPeaksCheck("Sicer",WkgDir,Config)){   
   
   WindowSize <- getSicerWindowSize(WkgDir,ConfigDirectory="Config")
   GapSize <- getSicerGapSize(WkgDir,ConfigDirectory="Config")      
   GenomeForSicer <- GetGenomeFromConfig(WkgDir,ConfigDirectory="Config")
   GenomeForSicer[tolower(GenomeForSicer) == "grch37"] <- "hg19" 
   GenomeForSicer[tolower(GenomeForSicer) == "hg18"] <- "hg18"
   GenomeForSicer[GenomeForSicer == "mm9"] <- "mm9"
  #"hg18","GRCh37","hg19","MM9","mm9"


   Pipeline <- getPipelinesPath("sicerpeakcallpipeline")
   workflowExec <- getWorkflowParam("Executable")
   javaExec <- getExecPath("java")
   mode <- tolower(getWorkflowParam("Mode"))
   PipelineBase <- GetPipelinebase()
   pythonExec <- getExecPath("python")    
   pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")
   rexec <- getExecPath("rexec")
   sicerexec <- getSicerExec()
   pythonlibs <- getLibraryPath("pythonlibs")

   
   ShiftSizeDefault <- 100
   
   dir.create(file.path(PLs@SicerDir,""),recursive=T,showWarnings=F)
    
   #genomeChrLengths <- GetChrLengthsFromConfig(WkgDir,ConfigDirectory="Config")
   SampleSheet <- SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe",]
   SamplesAndInputs <- matrix(data=c(gsub("\\.bam","",SampleSheet[,"Processed_bamFileName"]),gsub("\\.bam","",SampleSheet[match(SampleSheet[,"InputToUse"],SampleSheet[,"SampleName"],incomparable=NA),"Processed_bamFileName"]),as.vector(SampleSheet[,"CC_Fragment_Length"])),ncol=3,byrow=F,dimnames=list(NULL,c("Samples","Inputs","CC_Fraglength")))
   SamplesAndInputs <- SamplesAndInputs[!is.na(SamplesAndInputs[,"Inputs"]) & ! as.vector(SamplesAndInputs[,"Samples"]) %in% "No_Processed_Bam",]

   SamplesAndInputs2 <- matrix(data=c(gsub("\\.bam","",SampleSheet[,"Processed_bamFileName"]),gsub("\\.bam","",SampleSheet[match(SampleSheet[,"InputToUse"],SampleSheet[,"GenomicsID"],incomparable=NA),"Processed_bamFileName"]),as.vector(SampleSheet[,"CC_Fragment_Length"])),ncol=3,byrow=F,dimnames=list(NULL,c("Samples","Inputs","CC_Fraglength")))  
   SamplesAndInputs2 <- SamplesAndInputs2[!is.na(SamplesAndInputs2[,"Inputs"]) & ! as.vector(SamplesAndInputs2[,"Samples"]) %in% "No_Processed_Bam",]

   MakeUniquePeakCalls <- rbind(SamplesAndInputs,SamplesAndInputs2)
   MakeUniquePeakCalls <- MakeUniquePeakCalls[match(unique(MakeUniquePeakCalls[,1]),MakeUniquePeakCalls[,1]),]   
   SamplesAndInputs <- matrix(MakeUniquePeakCalls,ncol=3,byrow=F)

   SamplesAndInputs[SamplesAndInputs[,3] %in% "Too_few_Reads_To_Calculate" | SamplesAndInputs[,3] %in% "No_Information_Available" | is.na(SamplesAndInputs[,3]),3] <-  ShiftSizeDefault

  Variables  <- c(file.path(PLs@SicerDir,""),file.path(PLs@BamDir,""),WkgDir,GenomeForSicer,sicerexec,WindowSize,GapSize,0.01,pythonlibs,PipelineBase,pythonExec)
  names(Variables) <- c("SicerDirectory","BamDirectory","WorkingDirectory","GS","sicerexec","Window","GapSize","FDR","pythonlibs","pipelineBase","pythonExec")
  

  if(nrow(SamplesAndInputs) > 0){
    Specialisations  <- vector("list",length=nrow(SamplesAndInputs))
    for(i in 1:length(Specialisations)){
      Specialisations[[i]] <- c(as.vector(SamplesAndInputs[i,1]),as.vector(SamplesAndInputs[i,2]),as.vector(SamplesAndInputs[i,3]))
      names(Specialisations[[i]]) <- c("Test","Control","SS")  
    } 
    names(Specialisations) <-  paste(as.vector(SamplesAndInputs[,1]),"_SicerPeakCall",sep="") 
  
    #Pipeline = "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/SicerPeakCallingPipeline.xml"  
    PipeName <- "SS_SicerPeak"
    SSSicerPeak_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=75)
    saveXML(SSSicerPeak_WfMeta,file=file.path(PLs@WorkFlowDir,"SS_SicerPeak.xml"))
    system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,"SS_SicerPeak.xml"),sep=""),wait=T)
  }

 #######################
#  MacsFiles <- dir(path=PLs@MacsDir,pattern="*peaks.xls$",full.name=F)
#  MacsFilesFull <- dir(path=PLs@MacsDir,pattern="*peaks.xls$",full.name=T)  
#  PositivePeaks <- MacsFiles[-grep("negative_peaks.xls",MacsFiles)]
#  PositivePeaksFull <- MacsFilesFull[-grep("negative_peaks.xls",MacsFilesFull)]

  SampleSheet <- ReadAndLock(file.path(WkgDir,"SampleSheet.csv"),WkdDir,SAF=F,napTime=5)
  
  tryCatch({ 
    SampleSheetTemp <- SampleSheet


  SicerPeaksFull <-  dir(path=PLs@SicerDir,recursive=T,full.names=T)
    
  for(i in 1:nrow(SampleSheetTemp)){
	 if(any(file.path(PLs@SicerDir,gsub(".bam","",SampleSheetTemp[i,"Processed_bamFileName"]),paste(gsub(".bam","",SampleSheetTemp[i,"Processed_bamFileName"]),"-W",WindowSize,"-G",GapSize,"-islands-summary-FDR0.01",sep="")) %in% SicerPeaksFull)){
	   	SicerToRead <- SicerPeaksFull[SicerPeaksFull %in% file.path(PLs@SicerDir,gsub(".bam","",SampleSheetTemp[i,"Processed_bamFileName"]),paste(gsub(".bam","",SampleSheetTemp[i,"Processed_bamFileName"]),"-W",WindowSize,"-G",GapSize,"-islands-summary-FDR0.01",sep=""))]
 	    print(SicerToRead)
 	    if(file.info(SicerToRead)$size > 0 & file.exists(SicerToRead)){
		  DataIn <- read.delim(SicerToRead,sep="\t",comment.char="#")
		  SampleSheetTemp[i,"SicerPeaks"] <- nrow(DataIn)		
		  SampleSheetTemp[i,"Sicer_name"] <- SicerToRead
		  print(nrow(DataIn))
		  }
    }
  }
   SampleSheet <- SampleSheetTemp 	
  },warning = function(war){
    write.table(war,file.path(PLs@WorkFlowDir,paste(PipeName,"_WarningFile.txt",sep="")),sep="\t")
  
  },error = function(err){
    write.table(err,file.path(PLs@WorkFlowDir,paste(PipeName,"_ErrorFile.txt",sep="")),sep="\t")
  
  },finally = {
    WriteAndUnlock(SampleSheet,file.path(WkgDir,"SampleSheet.csv"))
  }
  
  )
}        
}
###########################




RunTPICsPeakCallerPipeline <- function(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipelineLocations,Config="Config"){
if(CallPeaksCheck("TPICs",WkgDir,Config)){      
   
   min_size <- getTPICSminSize(WkgDir,ConfigDirectory="Config")
   WideScale <- getTPICSwideSize(WkgDir,ConfigDirectory="Config")
   GenomeForTPICS <- GetGenomeFromConfig(WkgDir,ConfigDirectory="Config")
#   GenomeForTPICS[GenomeForTPICS == "grch37"] <- "hg19" 
#   GenomeForSicer[GenomeForSicer == "grch37"] <- "mm9"



   Pipeline <- getPipelinesPath("tpicspeakcallpipeline")
   workflowExec <- getWorkflowParam("Executable")
   javaExec <- getExecPath("java")
   mode <- tolower(getWorkflowParam("Mode"))
   PipelineBase <- GetPipelinebase()
   perl <- getExecPath("perl")    
   pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")
   rexec <- getExecPath("rexec")
   tpics <- getTPICsCustomExec()
   tpicszeta <- getTPICsZetaCustomExec() 
   tpicscoverage <- getTPICsCreateCoverageCustomExec()   
   bedtools <- getExecPath("bedtools")   
   rlibs <- getLibraryPath("rlibs")
   genomeChrLengths <- GetChrLengthsFromConfig(WkgDir,ConfigDirectory="Config")
       
   
   ShiftSizeDefault <- 100
   
   dir.create(file.path(PLs@TPICsDir,""),recursive=T,showWarnings=F)
    
   #genomeChrLengths <- GetChrLengthsFromConfig(WkgDir,ConfigDirectory="Config")
   SampleSheet <- SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe",]
   SamplesAndInputs <- matrix(data=c(gsub("\\.bam","",SampleSheet[,"Processed_bamFileName"]),gsub("\\.bam","",SampleSheet[match(SampleSheet[,"InputToUse"],SampleSheet[,"SampleName"],incomparable=NA),"Processed_bamFileName"]),as.vector(SampleSheet[,"CC_Fragment_Length"])),ncol=3,byrow=F,dimnames=list(NULL,c("Samples","Inputs","CC_Fraglength")))
   SamplesAndInputs <- SamplesAndInputs[!is.na(SamplesAndInputs[,"Inputs"]) & ! as.vector(SamplesAndInputs[,"Samples"]) %in% "No_Processed_Bam",]

   SamplesAndInputs2 <- matrix(data=c(gsub("\\.bam","",SampleSheet[,"Processed_bamFileName"]),gsub("\\.bam","",SampleSheet[match(SampleSheet[,"InputToUse"],SampleSheet[,"GenomicsID"],incomparable=NA),"Processed_bamFileName"]),as.vector(SampleSheet[,"CC_Fragment_Length"])),ncol=3,byrow=F,dimnames=list(NULL,c("Samples","Inputs","CC_Fraglength")))  
   SamplesAndInputs2 <- SamplesAndInputs2[!is.na(SamplesAndInputs2[,"Inputs"]) & ! as.vector(SamplesAndInputs2[,"Samples"]) %in% "No_Processed_Bam",]

   MakeUniquePeakCalls <- rbind(SamplesAndInputs,SamplesAndInputs2)
   MakeUniquePeakCalls <- MakeUniquePeakCalls[match(unique(MakeUniquePeakCalls[,1]),MakeUniquePeakCalls[,1]),]   
   SamplesAndInputs <- matrix(MakeUniquePeakCalls,ncol=3,byrow=F)

   SamplesAndInputs[SamplesAndInputs[,3] %in% "Too_few_Reads_To_Calculate" | SamplesAndInputs[,3] %in% "No_Information_Available" | is.na(SamplesAndInputs[,3]),3] <-  ShiftSizeDefault

  Variables  <- c(file.path(PLs@TPICsDir,""),file.path(PLs@BamDir,""),WkgDir,GenomeForTPICS,min_size,WideScale,bedtools,tpics,tpicszeta,tpicscoverage,rexec,perl,rlibs,genomeChrLengths,PipelineBase)
  names(Variables) <- c("TPICsDirectory","BamDirectory","WorkingDirectory","GT","min_size","WideScale","bedtools","tpics","tpicszeta","tpicscoverage","rexec","perl","rlibs","genomeChrLengths","pipelineBase")
  

  if(nrow(SamplesAndInputs) > 0){
    Specialisations  <- vector("list",length=nrow(SamplesAndInputs))
    for(i in 1:length(Specialisations)){
      Specialisations[[i]] <- c(as.vector(SamplesAndInputs[i,1]),as.vector(SamplesAndInputs[i,2]),as.vector(SamplesAndInputs[i,3]))
      names(Specialisations[[i]]) <- c("Test","Control","SS")  
    } 
    names(Specialisations) <-  paste(as.vector(SamplesAndInputs[,1]),"_TPICsPeakCall",sep="") 
  
    #Pipeline = "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/TPICSPeakCallingPipeline.xml"  
    PipeName <- "SS_TPICsPeak"
    SSTPICsPeak_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=75)
    saveXML(SSTPICsPeak_WfMeta,file=file.path(PLs@WorkFlowDir,"SS_TPICsPeak.xml"))
    system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,"SS_TPICsPeak.xml"),sep=""),wait=T)
  }


  SampleSheet <- ReadAndLock(file.path(WkgDir,"SampleSheet.csv"),WkdDir,SAF=F,napTime=5)

  tryCatch({ 
    SampleSheetTemp <- SampleSheet


  TPICsPeaksFull <-  dir(path=PLs@TPICsDir,recursive=T,full.names=T)
    
  for(i in 1:nrow(SampleSheetTemp)){
	 if(any(file.path(PLs@TPICsDir,gsub(".bam","",SampleSheetTemp[i,"Processed_bamFileName"]),paste(gsub(".bam","",SampleSheetTemp[i,"Processed_bamFileName"]),"_TPICS_Peaks.bed",sep="")) %in% TPICsPeaksFull)){
	   	TPICsToRead <- TPICsPeaksFull[TPICsPeaksFull %in% file.path(PLs@TPICsDir,gsub(".bam","",SampleSheetTemp[i,"Processed_bamFileName"]),paste(gsub(".bam","",SampleSheetTemp[i,"Processed_bamFileName"]),"_TPICS_Peaks.bed",sep=""))]
 	    print(TPICsToRead)
 	    if(file.info(TPICsToRead)$size > 0 & file.exists(TPICsToRead)){
		    DataIn <- read.delim(TPICsToRead,sep="\t",comment.char="#")
		    SampleSheetTemp[i,"TPICsPeaks"] <- nrow(DataIn)		
		    SampleSheetTemp[i,"TPICs_name"] <- TPICsToRead
		  }else{
		    SampleSheetTemp[i,"TPICsPeaks"] <- 0
        DataIn <- 0	
        SampleSheetTemp[i,"TPICs_name"] <- NA	
		    #SampleSheetTemp[i,"TPICs_name"] <- TPICsToRead		  		  
		  }
		  
		  print(nrow(DataIn))
    }
  }
    
   SampleSheet <- SampleSheetTemp 	
  },warning = function(war){
    write.table(war,file.path(PLs@WorkFlowDir,paste(PipeName,"_WarningFile.txt",sep="")),sep="\t")
  
  },error = function(err){
    write.table(err,file.path(PLs@WorkFlowDir,paste(PipeName,"_ErrorFile.txt",sep="")),sep="\t")
  
  },finally = {
    WriteAndUnlock(SampleSheet,file.path(WkgDir,"SampleSheet.csv"))
  }
  
  )
        return(SampleSheet)
}        
}
###########################


RunPeakProfilingPipeline <- function(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipelineLocations,Config="Config",Caller="Macs"){
if(CallProfilesCheck(Caller,WkgDir,Config)){
  if(Caller == "Macs"){
  PeakDirectory <- PLs@MacsDir
  PeakFileNameColumn <- colnames(SampleSheet) %in% "Macs_name"
  } 
  if(Caller == "Sicer"){
  PeakDirectory <- PLs@SicerDir
  PeakFileNameColumn <- colnames(SampleSheet) %in% "Sicer_name"  
  } 
  if(Caller == "TPICs"){
  PeakDirectory <- PLs@TPICsDir
  PeakFileNameColumn <- colnames(SampleSheet) %in% "TPICs_name"  
  }   
  BamDir <- PLs@BamDir
  

  Pipeline <- getPipelinesPath("peakprofilepipeline")
  workflowExec <- getWorkflowParam("Executable")
  javaExec <- getExecPath("java")
  mode <- tolower(getWorkflowParam("Mode"))
  PipelineBase <- GetPipelinebase()
  pythonExec <- getExecPath("python")    
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")
  rexec <- getExecPath("rexec")
  bedtools <- getExecPath("bedtools")
  rlibs <- getLibraryPath("rlibs")
  

  GenePos <- GetGenePosFromConfig(WkgDir,ConfigDirectory="Config")
  GeneSets <- GetGeneSetsFromConfig(WkgDir,ConfigDirectory="Config")
  
  
  genome <- GetGenomeFromConfig(WkgDir,ConfigDirectory="Config") 
  genomeChrLengths <- GetChrLengthsFromConfig(WkgDir,ConfigDirectory="Config")
  fastaFromGenome <- GetGenomeBuildFromConfig(WkgDir,ConfigDirectory="Config")
  Variables  <- c(WkgDir,file.path(PeakDirectory,""),file.path(BamDir,""),genome,genomeChrLengths,fastaFromGenome,bedtools,rexec,PipelineBase,rlibs,GenePos,GeneSets,Caller)
  names(Variables) <- c("WorkingDirectory","PeakDirectory","BamDirectory","genome","genomeFile","fastafile","bedtools","rexec","pipelineBase","rlibs","genepos","genesets","caller")
  
  JustOfInterest <- SampleSheet[!is.na(SampleSheet[,PeakFileNameColumn]) & !SampleSheet[,PeakFileNameColumn] %in% "NA",]

  dir.create(file.path(PeakDirectory,"PeakProfiles"),recursive=T,showWarnings=F)
  if(nrow(JustOfInterest) > 0){
    Specialisations  <- vector("list",length=nrow(JustOfInterest))
    for(i in 1:length(Specialisations)){
      Specialisations[[i]] <- c(gsub(".xls",".bed",as.vector(JustOfInterest[i,PeakFileNameColumn])),as.vector(JustOfInterest[i,"GenomicsID"]),as.vector(JustOfInterest[i,"Processed_bamFileName"]),36,as.vector(JustOfInterest[i,"GenomicsID"]))
      names(Specialisations[[i]]) <- c("peakbed","mergedpeakbed","bam","ReadLength","Test")  
    } 
    names(Specialisations) <-  paste(as.vector(JustOfInterest[,1]),paste("_",Caller,"PeakCall",sep=""),sep="") 
  
    #Pipeline = "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/PeakProfilingPipeline2.xml"  
    PipeName <- paste("SS_",Caller,"PeakProfile",sep="")
    SSPeakProfile_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=75)
    saveXML(SSPeakProfile_WfMeta,file=file.path(PLs@WorkFlowDir,paste("SS_",Caller,"PeakProfile.xml",sep="")))
    system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,paste("SS_",Caller,"PeakProfile.xml",sep="")),sep=""),wait=T)
  }
#if(Caller %in% "Macs"){
      SampleSheet <- ReadAndLock(file.path(WkgDir,"SampleSheet.csv"),WkdDir,SAF=F,napTime=5)
      
      tryCatch({ 
      SampleSheetTemp <- SampleSheet
        
      
      
      PeakGenes <- gsub("_peaks_Annotated.xls","",dir(path=file.path(PeakDirectory),pattern="*Annotated.xls$"))
      PeakGenesFullPath <- dir(path=file.path(PeakDirectory),pattern="*Annotated.xls$",full.names=T)
    for(i in 1:nrow(SampleSheetTemp)){
    	if(length(grep(gsub(".bwa.homo_sapiens.bam","",SampleSheetTemp[i,"Source_File"]),PeakGenes)) > 0){
    		MacsToAnnotate <- PeakGenesFullPath[grep(gsub(".bwa.homo_sapiens.bam","",SampleSheetTemp[i,"Source_File"]),PeakGenes)]
    		SampleSheetTemp[i,paste(Caller,"_Genes_Annotation",sep="")] <- MacsToAnnotate
    		print("Hello")
    	}
    	if(length(grep(SampleSheetTemp[i,"GenomicsID"],PeakGenes)) > 0){
        MacsToAnnotate <- PeakGenesFullPath[PeakGenes %in% SampleSheetTemp[i,"GenomicsID"]]
        SampleSheetTemp[i,paste(Caller,"_Genes_Annotation",sep="")] <- MacsToAnnotate
    		print("Hello")
    	}
    }
    
        PeakGenes <- gsub("_peaks_Annotated.summary","",dir(path=file.path(PeakDirectory),pattern="*Annotated.summary$"))
        PeakGenesFullPath <- dir(path=file.path(PeakDirectory),pattern="*Annotated.summary$",full.names=T)        
    for(i in 1:nrow(SampleSheetTemp)){
    	if(length(grep(gsub(".bwa.homo_sapiens.bam","",SampleSheetTemp[i,"Source_File"]),PeakGenes)) > 0){
    		MacsToAnnotate <- PeakGenesFullPath[grep(gsub(".bwa.homo_sapiens.bam","",SampleSheetTemp[i,"Source_File"]),PeakGenes)]
    		tryCatch({ 
    		DataInAnno <- read.delim(MacsToAnnotate,sep="\t",comment.char="#",header=T)[,1]
    		SampleSheetTemp[i,paste(Caller,"_Genes_In_Peaks",sep="")] <- DataInAnno[1]
    		SampleSheetTemp[i,paste(Caller,"_Peaks_In_Genes",sep="")] <- DataInAnno[2]
    		},error = function(err){
    		SampleSheetTemp[i,paste(Caller,"_Genes_In_Peaks",sep="")] <- 0
    		SampleSheetTemp[i,paste(Caller,"_Peaks_In_Genes",sep="")] <- 0
        })
    	}
    	if(length(grep(SampleSheetTemp[i,"GenomicsID"],PeakGenes)) > 0){
  		  tryCatch({ 
    		MacsToAnnotate <- PeakGenesFullPath[PeakGenes %in% SampleSheetTemp[i,"GenomicsID"]]
    		DataInAnno <- read.delim(MacsToAnnotate,sep="\t",comment.char="#",header=T)[,1]
    		SampleSheetTemp[i,paste(Caller,"_Genes_In_Peaks",sep="")] <- DataInAnno[1]
    		SampleSheetTemp[i,paste(Caller,"_Peaks_In_Genes",sep="")] <- DataInAnno[2]   		
    		},error = function(err){
    		SampleSheetTemp[i,paste(Caller,"_Genes_In_Peaks",sep="")] <- 0
    		SampleSheetTemp[i,paste(Caller,"_Peaks_In_Genes",sep="")] <- 0    		
    		})
    	}
    	
    }
    PeaksGO <- gsub("_GO_Results.txt","",dir(path=file.path(PeakDirectory),pattern="*_GO_Results.txt$"))
    PeaksGOFullPath <- dir(path=file.path(PeakDirectory),pattern="*_GO_Results.txt$",full.names=T)    
    for(i in 1:nrow(SampleSheetTemp)){
      if(length(grep(gsub(".bwa.homo_sapiens.bam","",SampleSheetTemp[i,"Source_File"]),PeaksGO)) > 0){
    		PeaksGOToAnnotateSum <- PeaksGOFullPath[PeaksGO %in% gsub(".bwa.homo_sapiens.bam","",SampleSheetTemp[i,"Source_File"])]
    		tryCatch({     		
    		DataInAnno <- read.delim(PeaksGOToAnnotateSum,sep="\t",comment.char="#",header=T)
    		SampleSheetTemp[i,paste(Caller,"_GO_Terms",sep="")] <- nrow(DataInAnno[DataInAnno[,"Adjusted_P_Value"] < 0.01 & !is.na(DataInAnno[,"Adjusted_P_Value"]),])
    		SampleSheetTemp[i,paste(Caller,"_GO_Table",sep="")] <- PeaksGOToAnnotateSum
    		},error = function(err){
    		SampleSheetTemp[i,paste(Caller,"_GO_Terms",sep="")] <- 0
    		SampleSheetTemp[i,paste(Caller,"_GO_Table",sep="")] <- NA
    		})
    		
    	}
    	if(length(grep(SampleSheetTemp[i,"GenomicsID"],PeaksGO)) > 0){
    		PeaksGOToAnnotateSum <- PeaksGOFullPath[PeaksGO %in% SampleSheetTemp[i,"GenomicsID"]]
    		tryCatch({ 
    		DataInAnno <- read.delim(PeaksGOToAnnotateSum,sep="\t",comment.char="#",header=T)
    		SampleSheetTemp[i,paste(Caller,"_GO_Terms",sep="")] <- nrow(DataInAnno[DataInAnno[,"Adjusted_P_Value"] < 0.01 & !is.na(DataInAnno[,"Adjusted_P_Value"]),])
    		SampleSheetTemp[i,paste(Caller,"_GO_Table",sep="")] <- PeaksGOToAnnotateSum
        },error = function(err){
    		SampleSheetTemp[i,paste(Caller,"_GO_Terms",sep="")] <- 0
    		SampleSheetTemp[i,paste(Caller,"_GO_Table",sep="")] <- NA
    		})
        
      }    	
    
    }
   SampleSheet <- SampleSheetTemp 	
  },warning = function(war){
    write.table(war,file.path(PLs@WorkFlowDir,paste(PipeName,"_WarningFile.txt",sep="")),sep="\t")
  
  },error = function(err){
    write.table(err,file.path(PLs@WorkFlowDir,paste(PipeName,"_ErrorFile.txt",sep="")),sep="\t")
  
  },finally = {
    WriteAndUnlock(SampleSheet,file.path(WkgDir,"SampleSheet.csv"))
  }
  
  )
    
#}


}
}

RunPeakMotifingPipeline <- function(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipelineLocations,Config="Config",Caller="Macs",NPeaks=500){
if(CallMotifsCheck(Caller,WkgDir,Config)){   
  if(Caller == "Macs"){
  PeakDirectory <- PLs@MacsDir
  PeakFileNameColumn <- colnames(SampleSheet) %in% "Macs_name"
  ColumnIndexForRank <- 5
  RankOrder <- "Rank"
  } 
  if(Caller == "Sicer"){
  PeakDirectory <- PLs@SicerDir
  PeakFileNameColumn <- colnames(SampleSheet) %in% "Sicer_name" 
  ColumnIndexForRank <- 7  
  RankOrder <- "RevRank"
  } 
  if(Caller == "TPICs"){
  PeakDirectory <- PLs@TPICsDir
  PeakFileNameColumn <- colnames(SampleSheet) %in% "TPICs_name"  
  ColumnIndexForRank <- 6
  RankOrder <- "RevRank"    
  }   
  
  Pipeline <- getPipelinesPath("motifpipeline")
  workflowExec <- getWorkflowParam("Executable")
  javaExec <- getExecPath("java")
  mode <- tolower(getWorkflowParam("Mode"))
  PipelineBase <- GetPipelinebase()
  pythonExec <- getExecPath("python")    
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")
  rexec <- getExecPath("rexec")
  meme <- getExecPath("meme")
  ame <- getExecPath("ame")
  rlibs <- getLibraryPath("rlibs")
  perllibs <- getLibraryPath("perllibs")
  pythonlibs <- getLibraryPath("pythonlibs")
    
  genomeChrLengths <- GetChrLengthsFromConfig(WkgDir,ConfigDirectory="Config")
   
  
  
  BamDir <- PLs@BamDir
  MemeFormatdb <- GetTFDB()
  genome <- GetGenomeFromConfig(WkgDir,ConfigDirectory="Config") 
  genomeChrLengths <- GetChrLengthsFromConfig(WkgDir,ConfigDirectory="Config")
  fastaFromGenome <- GetGenomeBuildFromConfig(WkgDir,ConfigDirectory="Config")
  Variables  <- c(WkgDir,file.path(PeakDirectory,""),file.path(BamDir,""),genome,genomeChrLengths,fastaFromGenome,NPeaks,RankOrder,ColumnIndexForRank,MemeFormatdb,ame,meme,rexec,PipelineBase,pythonExec,rlibs,perllibs,pythonlibs,genomeChrLengths)
  names(Variables) <- c("WorkingDirectory","PeakDirectory","BamDirectory","genome","genomeFile","Fasta","NPeaks","RankOrder","RankColumn","MotifDatabaseLocation","ame","meme","rexec","pipelineBase","python","rlibs","perllibs","pythonlibs","chrlengths")
  
  JustOfInterest <- SampleSheet[!is.na(SampleSheet[,PeakFileNameColumn]) & !SampleSheet[,PeakFileNameColumn] %in% "NA",]
if(nrow(JustOfInterest) > 0){
  #  dir.create(file.path(PeakDirectory,"PeakProfiles"),recursive=T)
  if(Caller == "Macs"){  
    Specialisations  <- vector("list",length=nrow(JustOfInterest))
    for(i in 1:length(Specialisations)){
      Specialisations[[i]] <- c("",gsub(".xls",".bed",basename(as.vector(JustOfInterest[i,PeakFileNameColumn]))),as.vector(JustOfInterest[i,"GenomicsID"]))
      names(Specialisations[[i]]) <- c("SamplePeakDir","PeakFileName","Test")  
    } 
    names(Specialisations) <-  paste(as.vector(JustOfInterest[,1]),paste("_",Caller,"PeakCall",sep=""),sep="") 
  }
  if(Caller %in% c("Sicer","TPICs")){  
    Specialisations  <- vector("list",length=nrow(JustOfInterest))
    for(i in 1:length(Specialisations)){
      Specialisations[[i]] <- c(gsub(".bam","",as.vector(JustOfInterest[i,"Processed_bamFileName"])),gsub(".xls",".bed",basename(as.vector(JustOfInterest[i,PeakFileNameColumn]))),as.vector(JustOfInterest[i,"GenomicsID"]))
      names(Specialisations[[i]]) <- c("SamplePeakDir","PeakFileName","Test")  
    } 
    names(Specialisations) <-  paste(as.vector(JustOfInterest[,1]),paste("_",Caller,"PeakCall",sep=""),sep="") 
  }
  
  
#    Pipeline = "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/MotifProcessPipeline2.xml"  
    PipeName <- paste("SS_",Caller,"Motifs",sep="")
    SSPeakProfile_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=75)
    saveXML(SSPeakProfile_WfMeta,file=file.path(PLs@WorkFlowDir,paste("SS_",Caller,"MotifProcess.xml",sep="")))
    system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,paste("SS_",Caller,"MotifProcess.xml",sep="")),sep=""),wait=T)
}
  
  SampleSheet <- ReadAndLock(file.path(WkgDir,"SampleSheet.csv"),WkgDir,SAF=F,napTime=5)

  tryCatch({ 
  SampleSheetTemp <- SampleSheet



  for(i in 1:nrow(SampleSheet)){
  	 SampleToLookFor <- SampleSheet[i,"GenomicsID"]
  	 KnownMotifFile <- dir(path=file.path(PeakDirectory,"Motif",SampleToLookFor,"Known"),pattern="*me.txt$",full.names=T)
  	 DataIn <- 0
    	 if(all(c(length(KnownMotifFile) > 0,file.info(KnownMotifFile)$size > 0))){
    		SampleSheet[i,paste(Caller,"_","Known_Motif_File",sep="")] <- KnownMotifFile
    		print(KnownMotifFile)
    		try(
    		DataIn <- read.delim(KnownMotifFile,sep=" ",comment.char="#",skip=11,header=F)
    		,silent=T)
    		if(all(DataIn != 0) & !is.null(DataIn)){
      		try(
      		SampleSheet[i,paste(Caller,"_","Known_Significant_Motifs",sep="")]  <- nrow(DataIn)
      		,silent=T)
    		}else{
    		try(
    		SampleSheet[i,paste(Caller,"_","Known_Significant_Motifs",sep="")]  <- 0
      		,silent=T)
    		}
        		
      }
      DataIn <- 0
    }

    for(i in 1:nrow(SampleSheet)){
    	SampleToLookFor <- as.vector(SampleSheet[i,"GenomicsID"])
    	print(SampleToLookFor)
    	DremeMotifFile <- dir(path=file.path(PeakDirectory,"Motif",SampleToLookFor,"Denovo","dreme_out"),pattern="*.xml$",full.names=T)
    	if(all(c(length(DremeMotifFile) > 0,file.info(DremeMotifFile)$size > 0))){
    		doc = xmlTreeParse(DremeMotifFile, useInternal = TRUE)
    		top = xmlRoot(doc)
    		NOfMotifs <- length(names(top[["motifs"]]))
    		SampleSheet[i,paste(Caller,"_","Dreme_Motif_File",sep="")] <- dir(path=file.path(PeakDirectory,"Motif",SampleToLookFor,"Denovo","dreme_out"),pattern="*.html$",full.names=T)
    		SampleSheet[i,paste(Caller,"_","Dreme_Significant_Motifs",sep="")] <- NOfMotifs
    	}
    }


      for(i in 1:nrow(SampleSheet)){
    	
      SampleToLookFor <- as.vector(SampleSheet[i,"GenomicsID"])    	
    	DremeTomTomFile <- dir(path=file.path(PeakDirectory,"Motif",SampleToLookFor,"Denovo","dreme_tomtom_out"),pattern="*.txt$",full.names=T)
      print(i)    	
      if(all(c(length(DremeTomTomFile) > 0,file.info(DremeTomTomFile)$size > 0))){
    		TempTomTom <- read.delim(DremeTomTomFile,sep="\t",header=T)
    		if(nrow(TempTomTom) > 0){
    		PerSampleTomTom <- cbind(SampleToLookFor,TempTomTom)
    		if(exists(as.character(bquote(SignificantTomTomDreme)))){
    			SignificantTomTomDreme <- rbind(SignificantTomTomDreme,PerSampleTomTom)
    		}else{
    			SignificantTomTomDreme <- PerSampleTomTom 
    		}
    			SampleSheet[i,paste(Caller,"_","Dreme_Significant_Motifs_Matched_To_Known",sep="")] <- length(unique(PerSampleTomTom[,2]))
    		}else{
    			SampleSheet[i,paste(Caller,"_","Dreme_Significant_Motifs_Matched_To_Known",sep="")] <- 0
    		}
    	}
}

      for(i in 1:nrow(SampleSheet)){
        print(i)
      	SampleToLookFor <- SampleSheet[i,"GenomicsID"]
      	MemeMotifFile <- dir(path=file.path(PeakDirectory,"Motif",SampleToLookFor,"Denovo","meme_out"),pattern="*.xml$",full.names=T)
      	if(all(c(length(MemeMotifFile) > 0,file.info(MemeMotifFile)$size > 0))){
      		doc = xmlTreeParse(MemeMotifFile, useInternal = TRUE)
      		top = xmlRoot(doc)
      		NOfMotifs <- length(names(top[["motifs"]]))
      		SampleSheet[i,paste(Caller,"_","meme_Motif_File",sep="")] <- dir(path=file.path(PeakDirectory,"Motif",SampleToLookFor,"Denovo","meme_out"),pattern="*.html$",full.names=T)
      		SampleSheet[i,paste(Caller,"_","meme_Significant_Motifs",sep="")] <- NOfMotifs	
        	}
        }

      for(i in 1:nrow(SampleSheet)){
      print(i)
      	SampleToLookFor <- SampleSheet[i,"GenomicsID"]
       	memeTomTomFile <- dir(path=file.path(PeakDirectory,"Motif",SampleToLookFor,"Denovo","meme_tomtom_out"),pattern="*.txt$",full.names=T)
      	if(all(c(length(memeTomTomFile) > 0,file.info(memeTomTomFile)$size > 0))){
      		TempTomTom <- read.delim(memeTomTomFile,sep="\t",header=T)
      		if(nrow(TempTomTom) > 0){
      		PerSampleTomTom <- cbind(SampleToLookFor,TempTomTom)
      		if(exists(as.character(bquote(SignificantTomTommeme)))){
      			SignificantTomTommeme <- rbind(SignificantTomTommeme,PerSampleTomTom)
      		}else{
      			SignificantTomTommeme <- PerSampleTomTom 
      		}
      		SampleSheet[i,paste(Caller,"_","meme_Significant_Motifs_Matched_To_Known",sep="")] <- length(unique(PerSampleTomTom[,2]))
      		}else{
      		SampleSheet[i,paste(Caller,"_","meme_Significant_Motifs_Matched_To_Known",sep="")] <- 0
      		}
    	 }
      }

   #SampleSheet	
  },warning = function(war){
    write.table("Warning",file.path(PLs@WorkFlowDir,paste(PipeName,"_WarningFile.txt",sep="")),sep="\t")
    return("warning")
  },error = function(err){
    write.table("Error",file.path(PLs@WorkFlowDir,paste(PipeName,"_ErrorFile.txt",sep="")),sep="\t")
    return("error")  
  },finally = {
    WriteAndUnlock(SampleSheet,file.path(WkgDir,"SampleSheet.csv"))
  }
  
  )
        
    	
}    	
    	
}


RunBetweenPeaksPipeline <- function(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipelineLocations,Config="Config",Caller="Macs"){
if(CallBetweenPeaksCheck(Caller,WkgDir,Config)){ 
  if(Caller == "Macs"){
  PeakDirectory <- PLs@MacsDir
  PeakFileNameColumn <- colnames(SampleSheet) %in% "Macs_name"
  ColumnIndexForRank <- 5
  RankOrder <- "Rank"
  } 
  if(Caller == "Sicer"){
  PeakDirectory <- PLs@SicerDir
  PeakFileNameColumn <- colnames(SampleSheet) %in% "Sicer_name" 
  ColumnIndexForRank <- 6  
  RankOrder <- "RevRank"
  } 
  if(Caller == "TPICs"){
  PeakDirectory <- PLs@TPICsDir
  PeakFileNameColumn <- colnames(SampleSheet) %in% "TPICs_name"  
  ColumnIndexForRank <- 5
  RankOrder <- "RevRank"    
  }   
  BamDir <- PLs@BamDir
  
    
  Pipeline <- getPipelinesPath("betweenpeakspipeline")
  workflowExec <- getWorkflowParam("Executable")
  javaExec <- getExecPath("java")
  mode <- tolower(getWorkflowParam("Mode"))
  PipelineBase <- GetPipelinebase()
  pythonExec <- getExecPath("python")    
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")
  rexec <- getExecPath("rexec")
  rlibs <- getLibraryPath("rlibs")


  
  
  genomeChrLengths <- GetChrLengthsFromConfig(WkgDir,ConfigDirectory="Config")
  Variables  <- c(WkgDir,file.path(PeakDirectory,""),genomeChrLengths,PipelineBase,rexec,rlibs)
  names(Variables) <- c("WorkingDirectory","PeakDirectory","lengths","pipelineBase","rexec","rlibs")
  
  system(paste("mkdir -p ",file.path(PeakDirectory,"Between_Peaks",""),sep=""),wait=T)  
  
  

    
  JustOfInterest <- SampleSheet[!is.na(SampleSheet[,PeakFileNameColumn]) & !SampleSheet[,PeakFileNameColumn] %in% "NA",]
  if(nrow(JustOfInterest) > 0){
  AllPeaks <-   as.vector(JustOfInterest[,PeakFileNameColumn])
  AllPeaksNames <-   as.vector(JustOfInterest[,"GenomicsID"])
  speicalNames <-   vector("character")
  Specialisations  <- vector("list")

    Count=0
      for(i in 1:length(AllPeaks)){
          for(k in 1:length(AllPeaks)){
          Count=Count+1
            Specialisations[[Count]] <- c(
            gsub(".xls",".bed",as.vector(AllPeaks[i])),
            gsub(".xls",".bed",as.vector(AllPeaks[k])),
            as.vector(file.path(PeakDirectory,"Between_Peaks",paste(as.vector(AllPeaksNames[i]),"Vs",as.vector(AllPeaksNames[k]),sep="_")))
            )
            names(Specialisations[[Count]]) <- c("SampleName","SampleName2","OutName")
            speicalNames[Count] <-  paste(as.vector(AllPeaksNames[i]),"Vs",as.vector(AllPeaksNames[k]),sep="_")
      }
      } 
      names(Specialisations) <-  paste(speicalNames,"_BetweenPeaks",sep="") 
      
      #Pipeline = "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/BetweenPeaksPipeline2.xml"  
      PipeName <- paste("SS_",Caller,"Between_Peaks",sep="")
      SSBetweenPeaks_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=75)
      saveXML(SSBetweenPeaks_WfMeta,file=file.path(PLs@WorkFlowDir,paste("SS_",Caller,"BetweenPeaksProcess.xml",sep="")))
      system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,paste("SS_",Caller,"BetweenPeaksProcess.xml",sep="")),sep=""),wait=T)
    }
    }
}      	

#

RunAcrossPeaksPipeline <- function(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipelineLocations,Config="Config"){

  if(length(colnames(SampleSheet[,which(grepl("Peaks$",colnames(SampleSheet)) & !grepl("Genes",colnames(SampleSheet)))]) > 1)){

  Pipeline <- getPipelinesPath("acrosspeakspipeline")
  workflowExec <- getWorkflowParam("Executable")
  javaExec <- getExecPath("java")
  mode <- tolower(getWorkflowParam("Mode"))
  PipelineBase <- GetPipelinebase()
  pythonExec <- getExecPath("python")    
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")
  rexec <- getExecPath("rexec")
  rlibs <- getLibraryPath("rlibs")


  


  BamDir <- PLs@BamDir
  genomeChrLengths <- GetChrLengthsFromConfig(WkgDir,ConfigDirectory="Config")
  Variables  <- c(WkgDir,rexec,PipelineBase,rlibs)
  names(Variables) <- c("WorkingDirectory","rexec","pipelineBase","rlibs")
  
  system(paste("mkdir -p ",file.path(WkgDir,"Peaks","Across_Peaks",""),sep=""),wait=T)  
  GeneralPeakDir <- file.path(WkgDir,"Peaks","Across_Peaks","")
  
  PeakColumns <- SampleSheet[,c(1,which(grepl("Peaks$",colnames(SampleSheet)) & !grepl("Genes",colnames(SampleSheet)))+1)]
  JustPeaks <- PeakColumns[apply(PeakColumns[,-1],1,function(x)any(!is.na(x))),]
  
  PeakCaller <- gsub("_name","",colnames(SampleSheet[,which(grepl("Peaks$",colnames(SampleSheet)) & !grepl("Genes",colnames(SampleSheet)))+1]))
  PeakCallerCombinations <- combn(PeakCaller,2)

  PeakCallerCol <- colnames(SampleSheet[,which(grepl("Peaks$",colnames(SampleSheet)) & !grepl("Genes",colnames(SampleSheet)))+1])
  PeakCallerColCombinations <- combn(PeakCallerCol,2)
  
  
  
  Specialisations <- list()
  namesForSpecials <- vector()  
  count = 0
  
  for(i in 1:ncol(PeakCallerColCombinations)){
    for(k in 1:nrow(JustPeaks)){
          Peaks1 <- gsub(".xls",".bed",as.vector(JustPeaks[k,colnames(JustPeaks) %in% PeakCallerColCombinations[1,i]]))
          Peaks2 <- gsub(".xls",".bed",as.vector(JustPeaks[k,colnames(JustPeaks) %in% PeakCallerColCombinations[2,i]]))      
          if(!is.na(Peaks1) & !is.na(Peaks2)){
          count <- count+1
            OutFile2 <- file.path(GeneralPeakDir,paste(as.vector(JustPeaks[k,1]),"_",paste(PeakCallerCombinations[1,i],PeakCallerCombinations[2,i],sep="_Vs_"),"_Res.txt",sep=""))
            OutFile1 <- file.path(GeneralPeakDir,paste(as.vector(JustPeaks[k,1]),"_",paste(PeakCallerCombinations[1,i],PeakCallerCombinations[2,i],sep="_Vs_"),".bed",sep=""))
            namesForSpecials[count]   <- paste(as.vector(JustPeaks[k,1]),PeakCallerCombinations[1,i],PeakCallerCombinations[2,i],sep="_")            
            Specialisations[[count]] <- c(Peaks1,Peaks2,OutFile1,OutFile2)
            names(Specialisations[[count]])  <- c("Peaks1","Peaks2","Outfile1","Outfile2")            
          }
          
    }
  
  }
  if(count > 0){
  
  
  
  names(Specialisations) <- namesForSpecials
      
      
      #Pipeline = "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/AcrossPeaksPipeline2.xml"  
      PipeName <- paste("SS_Across_PeakCalls",sep="")
      SSBetweenPeaks_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=75)
      saveXML(SSBetweenPeaks_WfMeta,file=file.path(PLs@WorkFlowDir,paste("SS_AcrossPeakCallsProcess.xml",sep="")))
      system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,paste("SS_AcrossPeakCallsProcess.xml",sep="")),sep=""),wait=T)
  } 
      PeakCallerCombos <- paste(PeakCallerCombinations[1,],PeakCallerCombinations[2,],sep="_Vs_")
      AllAcrossPeaksRes <- dir(path=file.path(WkgDir,"Peaks","Across_Peaks",""),pattern="_Res.txt")
      SampleSheet <- ReadAndLock(file.path(WkgDir,"SampleSheet.csv"),WkgDir,SAF=F,napTime=5)
      
      
      tryCatch({ 
          SampleSheetTemp <- SampleSheet
        
        for(i in 1:nrow(SampleSheetTemp)){
          for(k in 1:length(PeakCallerCombos)){
            Temp <- AllAcrossPeaksRes[grep(PeakCallerCombos[k],AllAcrossPeaksRes)]
            Temp2 <- Temp[grep(SampleSheetTemp[i,1],Temp)]
            if(length(Temp2) > 0){
              if(!any(is.na(Temp2))){
               
               
                Results <- read.delim(file.path(WkgDir,"Peaks","Across_Peaks",Temp2),sep="\t",h=F)
                SampleSheetTemp[i,paste("Jaccard_Score_Of_BasePair_Overlap",PeakCallerCombos[k],sep="_")] <- Results[8,2]
                SampleSheetTemp[i,paste("Jaccard_Score_Of_Peaks_Overlap",PeakCallerCombos[k],sep="_")] <- Results[7,2]            
                SampleSheetTemp[i,paste("Percent_Of_Intersection_Peaks_In",strsplit(PeakCallerCombos[k],"_Vs_")[[1]][1],"for",PeakCallerCombos[k],sep="_")] <- Results[5,2]
                SampleSheetTemp[i,paste("Percent_Of_Intersection_Peaks_In",strsplit(PeakCallerCombos[k],"_Vs_")[[1]][2],"for",PeakCallerCombos[k],sep="_")] <- Results[6,2] 
                SampleSheetTemp[i,paste("CoOccurence_BedFile_For",PeakCallerCombos[k],sep="_")] <- gsub("_Res.txt",".bed",file.path(WkgDir,"Peaks","Across_Peaks",Temp2))              
              }
            }
            #SampleSheetTemp[i,PeakCallerCombos[k]] 
          }      
        }
       SampleSheet <- SampleSheetTemp 	
      },warning = function(war){
        write.table(war,file.path(PLs@WorkFlowDir,paste(PipeName,"_WarningFile.txt",sep="")),sep="\t")
      
      },error = function(err){
        write.table(err,file.path(PLs@WorkFlowDir,paste(PipeName,"_ErrorFile.txt",sep="")),sep="\t")
      
      },finally = {
        WriteAndUnlock(SampleSheet,file.path(WkgDir,"SampleSheet.csv"))
      }
      
      )
  
      return(SampleSheet)
  
  }
}      	


TopPeaksByRank <- function(InFile,OutFile,RankColumn,RankOrder,NPeaks){
  TempPeaks <-   read.delim(InFile,sep="\t",header=F)
  if(RankOrder %in% "RevRank"){
     NoOfPeaks <- nrow(TempPeaks)
     if(NoOfPeaks < NPeaks){NPeaks <- NoOfPeaks}
     ActualPeaks <- TempPeaks[order(TempPeaks[,RankColumn],decreasing=F),][1:NPeaks,]  
  }
  if(RankOrder %in% "Rank"){
     NoOfPeaks <- nrow(TempPeaks)
     if(NoOfPeaks < NPeaks){NPeaks <- NoOfPeaks}
     ActualPeaks <- TempPeaks[order(TempPeaks[,RankColumn],decreasing=T),][1:NPeaks,]          
  }
  write.table(ActualPeaks,OutFile,sep="\t",row.names=F,col.names=F,quote=F)              
}

RunMainPipeline <- function(WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipelineLocations,Config="Config"){
   
#  genome <- GetGenomeFromConfig(WkgDir,ConfigDirectory=Config)
  #Pipeline = "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/MainPipeline2.xml"
  Pipeline <- getPipelinesPath("mainpipeline")
  workflowExec <- getWorkflowParam("Executable")
  R_executable <- getExecPath("rexec")  
  javaExec <- getExecPath("java")
  mode <- tolower(getWorkflowParam("Mode"))
  PipelineBase <- GetPipelinebase()
  rlibs <- getLibraryPath("rlibs")
  
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")
  WorkFlowDirectory <- PLs@WorkFlowDir
  Config2 <- file.path(WkgDir,Config,"config.ini")
  #R_executable <- "/lustre/mib-cri/carrol09/Work/MyPipe/R/R-2.15.0/bin/Rscript"
  RandomTrackerNumber <-  getRandString()
  Variables  <- c(WorkFlowDirectory,RandomTrackerNumber,Config2,R_executable,PipelineBase,rlibs)
  names(Variables) <- c("WorkflowDir","PathwayTracker","Config","R_Executable","pipelineBase","rlibs")
  Specialisations <- vector("list")
  #Specialisations[[1]] <- "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/20121114_MontoyaAR_DN_Hes6ChIP/bamFiles/ID-LNCAP-LM-HES6-BICALUTAMIDE-INPUT-D26.bwa.bam"
  
  PipeName <- "Main"
  Main_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=75)
  saveXML(Main_WfMeta,file=file.path(PLs@WorkFlowDir,"Main.xml"))
  system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,"Main.xml"),sep=""),wait=T)

}

RunReportingPipeline <- function(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipelineLocations,Config="Config"){


  Pipeline <- getPipelinesPath("mainreportpipeline")
  workflowExec <- getWorkflowParam("Executable")
  R_executable <- getExecPath("rexec")  
  javaExec <- getExecPath("java")
  mode <- tolower(getWorkflowParam("Mode"))
  PipelineBase <- GetPipelinebase()
  rlibs <- getLibraryPath("rlibs")
  
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")


  speicalNames <-   vector("character")
  Specialisations  <- vector("list")
  Specialisations[[1]] <- "SampleSheet.csv"
  names(Specialisations[[1]]) <- c("SampleSheet")
  names(Specialisations) <-  "ReportMaker"
  Variables <- c(WkgDir,R_executable,PipelineBase,rlibs)
  names(Variables) <- c("WorkingDirectory","rexec","pipelineBase","rlibs")
      
  #Pipeline = "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/MainReportPipeline2.xml"  
  PipeName <- paste("SS_MakeReport",sep="")
  SSBetweenPeaks_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=75)
  saveXML(SSBetweenPeaks_WfMeta,file=file.path(PLs@WorkFlowDir,"SS_MakeReport.xml"))
  system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,paste("SS_MakeReport.xml",sep="")),sep=""),wait=T)

}      	



newXMLNode <- function (name, ..., attrs = NULL, namespace = character(), namespaceDefinitions = character(),
    doc = NULL, .children = list(...), parent = NULL, at = NA,
    cdata = FALSE, suppressNamespaceWarning = getOption("suppressXMLNamespaceWarning",
        FALSE), sibling = NULL, addFinalizer = NA, noNamespace = length(namespace) ==
        0 && !missing(namespace))
{
    if (length(attrs)) {
        ids = names(attrs)
        attrs = structure(as(attrs, "character"), names = ids)
        i = grep("^xmlns", names(attrs))
        if (length(i)) {
            warning("Don't specify namespace definitions via 'attrs'; use namespaceDefinitions")
            namespace = c(namespace, structure(attrs[i], names = gsub("^xmlns:",
                "", names(attrs)[i])))
            attrs = attrs[-i]
        }
    }
    else attrs = character()
    ns = character()
    name = strsplit(name, ";")[[1]]
    if (length(name) == 2) {
        ns = name[1]
        name = name[2]
        noNamespace = FALSE
    }
    if (is.list(parent)) {
        if (length(parent) < 1 || !(is(parent[[1]], "XMLInternalElementNode") ||
            is(parent[[1]], "XMLInternalDocument")))
            stop("incorrect value for parent")
        parent = parent[[1]]
    }
    if (missing(doc) && !missing(parent) && inherits(parent,
        "XMLInternalDocument")) {
        doc = parent
        parent = NULL
    }
    if (is.null(doc) && !is.null(parent)) {
        doc = if (inherits(parent, "XMLInternalDocument"))
            parent
        else .Call("R_getXMLNodeDocument", parent, PACKAGE = "XML")
    }
    node <- .Call("R_newXMLNode", as.character(name), character(),
        character(), doc, namespaceDefinitions, addFinalizer,
        PACKAGE = "XML")
    if (!is.null(sibling))
        addSibling(sibling, node, after = as.logical(at))
    else if (!is.null(parent))
        addChildren(parent, node, at = at)
    if (TRUE) {
        nsDefs = lapply(seq(along = namespaceDefinitions), function(i) newNamespace(node,
            namespaceDefinitions[[i]], names(namespaceDefinitions)[i],
            set = FALSE))
        if (length(namespaceDefinitions))
            names(nsDefs) = if (length(names(namespaceDefinitions)))
                names(namespaceDefinitions)
            else ""
    }
    else nsDefs = xmlNamespaceDefinitions(node)
    addAttributes(node, .attrs = attrs, suppressNamespaceWarning = suppressNamespaceWarning)
    if (is(namespace, "XMLNamespaceRef")) {
        setInternalNamespace(node, namespace)
    }
    else if (is.na(noNamespace) || !noNamespace) {
        ns =  XML:::getNodeNamespace(ns, nsDefs, node, namespace, noNamespace,
            namespaceDefinitions, parent, suppressNamespaceWarning)
        if (is.null(ns))
            .Call("R_setNamespaceFromAncestors", node, PACKAGE = "XML")
    }
    if (length(ns) && (inherits(ns, c("XMLNamespaceRef", "XMLNamespaceDeclaration")) ||
        (is.character(ns) && ns != "")))
        setXMLNamespace(node, ns)
    if (length(.children)) {
        if (!is.list(.children))
            .children = list(.children)
        addChildren(node, kids = .children, cdata = cdata, addFinalizer = addFinalizer)
    }
    node
}




RunDownSamplePipeline <- function(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipelineLocations,Config="Config"){
if(DownSampleCheck(WkgDir,Config)){
   

  Pipeline <- getPipelinesPath("downsamplepipeline")
  workflowExec <- getWorkflowParam("Executable")
  javaExec <- getExecPath("java")
  mode <- tolower(getWorkflowParam("Mode"))
  PipelineBase <- GetPipelinebase()   
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")
  picardExec <- getExecPath("picard")    

   
    if(GetDupFlag(WkgDir,Config) == "True"){
      TotalColumn <-  "Unique"
    }else{
      TotalColumn <-  "Filtered"   
    }
   
    #genomeChrLengths <- GetChrLengthsFromConfig(WkgDir,ConfigDirectory="Config")
   SampleSheet <- SampleSheet[SampleSheet[,"Analysis_State"] %in% "RunMe",]
   #SamplesAndInputs <- matrix(data=c(gsub("\\.bam","",SampleSheet[,"Processed_bamFileName"]),gsub("\\.bam","",SampleSheet[match(SampleSheet[,"InputToUse"],SampleSheet[,"SampleName"],incomparable=NA),"Processed_bamFileName"]),as.vector(SampleSheet[,TotalColumn]),as.vector(SampleSheet[match(SampleSheet[,"InputToUse"],SampleSheet[,"SampleName"],incomparable=NA),TotalColumn])),ncol=4,byrow=F,dimnames=list(NULL,c("Samples","Inputs","TotalSample","TotalInput")))
   #SamplesAndInputs <- SamplesAndInputs[!is.na(SamplesAndInputs[,"Inputs"]) & ! as.vector(SamplesAndInputs[,"Samples"]) %in% "No_Processed_Bam",]

   SamplesAndInputs2 <- matrix(data=c(gsub("\\.bam","",SampleSheet[,"Processed_bamFileName"]),gsub("\\.bam","",SampleSheet[match(SampleSheet[,"InputToUse"],SampleSheet[,"GenomicsID"],incomparable=NA),"Processed_bamFileName"]),as.vector(SampleSheet[,TotalColumn]),as.vector(SampleSheet[match(SampleSheet[,"InputToUse"],SampleSheet[,"GenomicsID"],incomparable=NA),TotalColumn])),ncol=4,byrow=F,dimnames=list(NULL,c("Samples","Inputs","TotalSample","TotalInput")))
   SamplesAndInputs2 <- SamplesAndInputs2[!is.na(SamplesAndInputs2[,"Inputs"]) & ! as.vector(SamplesAndInputs2[,"Samples"]) %in% "No_Processed_Bam",]

   MakeUniquePeakCalls <- SamplesAndInputs2
   MakeUniquePeakCalls <- MakeUniquePeakCalls[match(unique(MakeUniquePeakCalls[,1]),MakeUniquePeakCalls[,1]),]   
   SamplesAndInputs <- matrix(MakeUniquePeakCalls,ncol=4,byrow=F)
#   SamplesAndInputs2 <- SamplesAndInputs
#   SamplesAndInputs2[,3] <- SamplesAndInputs[,4]
#   SamplesAndInputs2[,4] <- SamplesAndInputs[,3]
#   SamplesAndInputs <- SamplesAndInputs2   
   
   Prob <- as.numeric(as.vector(SamplesAndInputs[,3]))/as.numeric(as.vector(SamplesAndInputs[,4]))
#   SamplesAndInputs[SamplesAndInputs[,3] %in% "Too_few_Reads_To_Calculate" | SamplesAndInputs[,3] %in% "No_Information_Available" |  is.na(SamplesAndInputs[,3]),3] <-  ShiftSizeDefault
   SamplesAndInputs <- cbind(SamplesAndInputs,Prob)
   SamplesAndInputs <- SamplesAndInputs[as.numeric(as.vector(SamplesAndInputs[,5])) < 1,]
   SamplesAndInputs <- SamplesAndInputs[!is.na(SamplesAndInputs[,1]),]
   Variables  <- c(file.path(PLs@BamDir,""),javaExec,picardExec)
   names(Variables) <- c("BamDirectory","java","picard")
  

  if(nrow(SamplesAndInputs) > 0){
    Specialisations  <- vector("list",length=nrow(SamplesAndInputs))
    for(i in 1:length(Specialisations)){
      Specialisations[[i]] <- c(as.vector(SamplesAndInputs[i,1]),as.vector(SamplesAndInputs[i,2]),as.vector(SamplesAndInputs[i,5]))
      names(Specialisations[[i]]) <- c("Sample","Input","prob")  
    } 
    names(Specialisations) <-  paste(as.vector(SamplesAndInputs[,1]),"_Downsample",sep="") 
  
    #Pipeline = "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/MacsPeakCallingPipeline.xml"  
    PipeName <- "SS_DownSample"
    SSMacsPeak_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=75)
    saveXML(SSMacsPeak_WfMeta,file=file.path(PLs@WorkFlowDir,"SS_Downsample.xml"))
    system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,"SS_Downsample.xml"),sep=""),wait=T)
}
  AllDownSampledFiles <- dir(PLs@BamDir,pattern="downsampled.*bam$",full.name=F)
  AllDownSampledFiles <-  AllDownSampledFiles[-grep("temp\\.bam",AllDownSampledFiles)]
 
  #SampleSheet <- ReadAndLock(file.path(WkgDir,"SampleSheet.csv"),WkdDir,SAF=F,napTime=5)

  for(i in 1:nrow(SampleSheet)){
	 if(gsub(".bam","",SampleSheet[i,"Processed_bamFileName"]) %in% gsub(".bam$","",gsub(".*downsampled_by_","",basename(AllDownSampledFiles)))){
	   	SampleSheet[i,"InputFileNameToUse"] <- AllDownSampledFiles[gsub(".bam$","",gsub(".*downsampled_by_","",basename(AllDownSampledFiles))) %in% gsub(".bam","",SampleSheet[i,"Processed_bamFileName"])]
		  
      }
	}
}
        return(SampleSheet)
}        


findSamples <- function(SampleSheet,SLXID=NULL,SampleName=NULL,Project=NULL){

SamplesFromProjects <- NULL
if(!is.null(Project)){
#Project <- c("20110623_ZecchiniV_DN_oxfChIP")
  for(i in 1:length(Project)){
    new <- getURL(paste("http://genomicsequencing.cruk.cam.ac.uk:8080/glsintapi/projectSampleDetails?projectName=",Project[i],sep=""))
    Projects <- xmlToList(new)
    TempSN <- unlist(Projects[rownames(Projects) == "sample",])
    SamplesFromProjects <- c(SamplesFromProjects,TempSN[names(TempSN) == "sample.name"])
  }
}
  
SampleName <- c(SampleName,as.vector(SamplesFromProjects))
library(XML)
library(RCurl)
library(RJSONIO)
SLXID <- c("SLX-7483")

LimsLocations <- NULL
if(!is.null(SLXID)){
  for(i in 1:length(SLXID)){
    new <- getURL(paste("http://genomicsequencing.cruk.cam.ac.uk:8080/glsintapi/runsContainingLibraries?slxId=",SLXID[i],sep=""))
    Runs <- xmlToList(new)
    slxidRuns <- vector("character")
    slxidTypes <- vector("character")
  
    TempCount <- which(Runs[1,] != "MiSeq Run")
  
    for(k in 1:length(TempCount)){
       slxidRuns[k] <- Runs[,TempCount[k]]$runFolder
       slxidTypes[k] <- Runs[,TempCount[k]]$runType
    }
    Temp <- cbind(rep(SLXID[i],length(slxidRuns)),slxidRuns,slxidTypes)
    #Temp <- Temp[-1,]
    colnames(Temp)  <- c("SLXID","Runs","Types")
  
     TempFull <- Temp
  
    for(j in 1:nrow(TempFull)){
       Files <- xmlToList(getURL(paste("http://genomicsequencing.cruk.cam.ac.uk:8080/glsintapi/fullDetailsOfRun?runId=",TempFull[j,"Runs"],sep="")))
       Lanes <- Files[rownames(Files) == "flowcell",][[1]][-c(1,length(Files[rownames(Files) == "flowcell",][[1]]))]
       incount = 0
       Temp2Full <- NULL
       for(l in 1:length(Lanes)){
         if(Lanes[[l]]$slxId == SLXID[i]){
            incount = incount+1
            print(l)
            JustSamples <- Lanes[[l]][names(Lanes[[l]]) == "sample"]
            vecSampleLocations <-  vector("character")
            vecSampleNames <-  vector("character")
            if(any(names(unlist(JustSamples)) == "sample.file.url")){
              for(f in 1:length(JustSamples)){
                JustSamples[[f]]$name
                vecSample <- matrix(unlist(JustSamples[[f]][names(JustSamples[[f]])== "file"]),nrow=length(unlist(JustSamples[[f]][names(JustSamples[[f]])== "file"]))/4,byrow=T)
                if(any(names(JustSamples[[f]]) == "file")){
                  if(any(grepl("Read 1 FASTQ",vecSample[,1]))){
                    vecSampleLocations[f] <- vecSample[grep("Read 1 FASTQ",vecSample[,1]),2]
                    vecSampleNames[f] <- JustSamples[[f]]$name
                  }
                }
  
  
              }
              Temp2 <- cbind(matrix(rep(TempFull[j,],length(vecSampleNames)),ncol=ncol(TempFull),byrow=T),rep(Lanes[[l]]$lane,length(vecSampleNames)),vecSampleNames,vecSampleLocations)
                #Temp <- Temp[-1,]
              colnames(Temp2)  <- c("SLXID","Run","Type","Lane","SampleNames","Location")
              if(incount > 1){
                 Temp2Full <- rbind(Temp2Full,Temp2)
              }else{
                 Temp2Full <- Temp2
              }
           }
         }
         
       }

          if(j > 1){
            Temp3 <- rbind(Temp3,Temp2Full)
          }else{
            Temp3 <- Temp2Full
          }

       
    }
        if(i > 1){
          LimsLocations <- rbind(LimsLocations,Temp3)
        }else{
          LimsLocations <- Temp3
        }
  
  
  }
}

SLXIDLimsLocations <- LimsLocations
LimsLocations <- NULL

if(!is.null(SampleName)){
  for(i in 1:length(SampleName)){
    new <- getURL(paste("http://genomicsequencing.cruk.cam.ac.uk:8080/glsintapi/runsContainingSamples?sampleName=",SampleName[i],sep=""))
    Runs <- xmlToList(new)
    samplenameRuns <- vector("character")
    samplenameTypes <- vector("character")
  
    TempCount <- which(Runs[1,] != "MiSeq Run")
  
    for(k in 1:length(TempCount)){
       samplenameRuns[k] <- Runs[,TempCount[k]]$runFolder
       samplenameTypes[k] <- Runs[,TempCount[k]]$runType
    }
    Temp <- cbind(rep("",length(samplenameRuns)),samplenameRuns,samplenameTypes)
    #Temp <- Temp[-1,]
    colnames(Temp)  <- c("SLXID","Runs","Types")
  
     TempFull <- Temp
  
    for(j in 1:nrow(TempFull)){
       Files <- xmlToList(getURL(paste("http://genomicsequencing.cruk.cam.ac.uk:8080/glsintapi/fullDetailsOfRun?runId=",TempFull[j,"Runs"],sep="")))
       Lanes <- Files[rownames(Files) == "flowcell",][[1]][-c(1,length(Files[rownames(Files) == "flowcell",][[1]]))]
       incount = 0
       Temp2Full <- NULL
       for(l in 1:length(Lanes)){
         #if(Lanes[[l]]$slxId == SLXID[i]){
            incount = incount+1
            print(l)
            
            LaneSLXID <- Lanes[[l]]$slxId
            
            JustSamples <- Lanes[[l]][names(Lanes[[l]]) == "sample"]
            vecSampleLocations <-  vector("character")
            vecSampleNames <-  vector("character")
            if(any(names(unlist(JustSamples)) == "sample.file.url")){
              for(f in 1:length(JustSamples)){
                if(JustSamples[[f]]$name == SampleName[i]){
                  vecSample <- matrix(unlist(JustSamples[[f]][names(JustSamples[[f]])== "file"]),nrow=4,byrow=T)
                  if(any(names(JustSamples[[f]]) == "file")){
                    if(any(grepl("Read 1 FASTQ",vecSample[,1]))){
                      vecSampleLocations[f] <- vecSample[grep("Read 1 FASTQ",vecSample[,1]),2]
                      
                      print("Hello")
                      vecSampleNames[f] <- JustSamples[[f]]$name
                      print(JustSamples[[f]]$name)
                      TempFull[j,1] <- LaneSLXID
                    }
                  }
               }
  
              }
              Temp2 <- cbind(matrix(rep(TempFull[j,],length(vecSampleNames)),ncol=ncol(TempFull),byrow=T),rep(Lanes[[l]]$lane,length(vecSampleNames)),vecSampleNames,vecSampleLocations)
                #Temp <- Temp[-1,]
              colnames(Temp2)  <- c("SLXID","Run","Type","Lane","SampleNames","Location")
              if(incount > 1){
                 Temp2Full <- rbind(Temp2Full,Temp2)
              }else{
                 Temp2Full <- Temp2
              }
           }
         #}
         
       }

          if(j > 1){
            Temp3 <- rbind(Temp3,Temp2Full)
          }else{
            Temp3 <- Temp2Full
          }

       
    }
        if(i > 1){
          LimsLocations <- rbind(LimsLocations,Temp3)
        }else{
          LimsLocations <- Temp3
        }
  
  
  }
}
  
SampleLimsLocations <- LimsLocations
AllLimsLocations <- rbind(SampleLimsLocations,SLXIDLimsLocations)
GenomicsID <- gsub(".fq.gz","",basename(AllLimsLocations[,"Location"]))
SampleName <- AllLimsLocations[,"SampleNames"]
Lane <- AllLimsLocations[,"Lane"]
Run <- AllLimsLocations[,"Run"]
FQLocations <- gsub("8080/solexa","",gsub("http://","",AllLimsLocations[,"Location"]))
SampleSheetTemp <- matrix(nrow=length(GenomicsID),ncol=ncol(SampleSheet))
colnames(SampleSheetTemp) <- colnames(SampleSheet)
SampleSheetTemp[,"GenomicsID"] <- GenomicsID
SampleSheetTemp[,"SampleName"] <- SampleName
SampleSheetTemp[,"Analysis_State"] <- "RunMe"
SampleSheetTemp[,"Lane"] <- Lane
SampleSheetTemp[,"Run"] <- Run
SampleSheetTemp[,"FQLocation"] <- FQLocations
SampleSheet <- rbind(SampleSheet,SampleSheetTemp)
if(nrow(SampleSheet) > 1){
  SampleSheet <- SampleSheet[!apply(SampleSheet,1,function(x)all(is.na(x))),]
}
return(SampleSheet)

}







RunConfigureAnnotationPipeline <- function(WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipeLineLocations,Config="Config"){


  Pipeline <- getPipelinesPath("configureannotationpipeline")
  workflowExec <- getWorkflowParam("Executable")
  javaExec <- getExecPath("java")
  mode <- tolower(getWorkflowParam("Mode"))
  PipelineBase <- GetPipelinebase()
  bwaExec <- getExecPath("bwa")
  pythonExec <- getExecPath("python")
  samtoolsExec <- getExecPath("samtools")  
  picardExec <- getExecPath("picard")    
  pythonlibs <- getLibraryPath("pythonlibs")

  GenomeBuild <- GetGenomeBuildFromConfig(WkgDir,ConfigDirectory=Config)

  
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")


  
  Variables <- c(WkgDir,bwaExec,pythonExec,javaExec,samtoolsExec,PipelineBase,picardExec,pythonlibs)
  names(Variables) <- c("WorkingDirectory","bwa","python","java","samtools","pipelineBase","picard","pythonlibs")

      Specialisations <- list(c(GenomeBuild))
      names(Specialisations[[length(Specialisations)]]) <- c("Genome")  
     names(Specialisations) <- seq(1,length(Specialisations))  
#     Pipeline <- "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/ReAlignPipeline.xml"
     PipeName <- "ConfigureAnnotation"
     SSRealignment_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=MaxJobs)
     saveXML(SSRealignment_WfMeta,file=file.path(PLs@WorkFlowDir,"ConfigureAnnotation.xml"))
     system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,"ConfigureAnnotation.xml"),sep=""),wait=T)
   
     Res <- readIniFile(file.path(Config,"config.ini"))
     ConfigOut <- file.path(Config,"config.ini")
     genome <- tolower(GetGenomeFromConfig(WkgDir,ConfigDirectory=Config))
     if(file.exists(paste(GenomeBuild,".sam",sep="")) & file.exists(paste(GenomeBuild,"_ChrLengths.txt",sep=""))){  
        SequenceDictionaryToUpdate <- Res[Res[,1] %in% "Sequence Dictionary" & Res[,2] %in% genome,3] <- paste(GenomeBuild,".sam",sep="")
        ChromosomeLengthsToUpdate <- Res[Res[,1] %in% "Chromosome Lengths" & Res[,2] %in% genome,3] <- paste(GenomeBuild,"_ChrLengths.txt",sep="")
     }
     
      AllSections <- unique(Res[,1])
      file.rename(ConfigOut,gsub(".ini","preupdatedanno.ini",ConfigOut))
      file.create(ConfigOut,showWarnings=T)
      for(i in 1:length(AllSections)){
        write.table(paste("[",AllSections[i],"]",sep=""),ConfigOut,col.names=F,row.names=F,quote=F,append=T)
        VariablesInSection <- as.vector(Res[Res[,1] %in% AllSections[i],2])
        for(j in 1:length(VariablesInSection)){
             Variable <- VariablesInSection[j]
             Value  <- as.vector(Res[Res[,1] %in% AllSections[i] & Res[,2] %in% VariablesInSection[j],3])
             write.table(paste(Variable,Value,sep=" = "),ConfigOut,col.names=F,row.names=F,quote=F,append=T)      
        }
        write.table("",ConfigOut,col.names=F,row.names=F,quote=F,append=T)      
      }

     
     
}



RunConfigureAnnotationPipeline2 <- function(WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipeLineLocations,Config="Config"){


  Pipeline <- getPipelinesPath("configureannotationpipeline")
  workflowExec <- getWorkflowParam("Executable")
  javaExec <- getExecPath("java")
  mode <- tolower(getWorkflowParam("Mode"))
  PipelineBase <- GetPipelinebase()
  bwaExec <- getExecPath("bwa")
  pythonExec <- getExecPath("python")
  samtoolsExec <- getExecPath("samtools")  
  picardExec <- getExecPath("picard")    
  gtftobedExec <- getExecPath("gtftobed")
  perlExec <- getExecPath("perl")
  rexec <- getExecPath("rexec")  
  
  pythonlibs <- getLibraryPath("pythonlibs")

  GenomeBuild <- GetGenomeBuildFromConfig(WkgDir,ConfigDirectory=Config)

  
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")

  GenePos <- GetGenePosFromConfig(WkgDir,ConfigDirectory="Config")
  #GeneSets <- GetGeneSetsFromConfig(WkgDir,ConfigDirectory="Config")
  if(length(grep(".gtf:",GenePos)) > 0){
     gtf <- strsplit(GenePos,":")$value[1]
     transcriptinfo <- strsplit(GenePos,":")$value[2]
     Pipeline <-  gsub(".xml","2.xml",Pipeline)
  }else
  {
     transcriptinfo <- ""
     gtf <- ""
  }
  
  Variables <- c(WkgDir,bwaExec,pythonExec,javaExec,samtoolsExec,PipelineBase,picardExec,pythonlibs,rexec,perlExec,gtftobedExec)
  names(Variables) <- c("WorkingDirectory","bwa","python","java","samtools","pipelineBase","picard","pythonlibs","rexec","perl","gtftobed")

      Specialisations <- list(c(GenomeBuild,gtf,transcriptinfo))
      names(Specialisations[[length(Specialisations)]]) <- c("Genome","gtf","transcriptinfo")  
     names(Specialisations) <- seq(1,length(Specialisations))  
#     Pipeline <- "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/ReAlignPipeline.xml"
     PipeName <- "ConfigureAnnotation"
     SSRealignment_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=MaxJobs)
     saveXML(SSRealignment_WfMeta,file=file.path(PLs@WorkFlowDir,"ConfigureAnnotation.xml"))
     system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,"ConfigureAnnotation.xml"),sep=""),wait=T)
   
     Res <- readIniFile(file.path(Config,"config.ini"))
     ConfigOut <- file.path(Config,"config.ini")
     genome <- tolower(GetGenomeFromConfig(WkgDir,ConfigDirectory=Config))
     if(file.exists(paste(GenomeBuild,".sam",sep="")) & file.exists(paste(GenomeBuild,"_ChrLengths.txt",sep=""))){  
        SequenceDictionaryToUpdate <- Res[Res[,1] %in% "Sequence Dictionary" & Res[,2] %in% genome,3] <- paste(GenomeBuild,".sam",sep="")
        ChromosomeLengthsToUpdate <- Res[Res[,1] %in% "Chromosome Lengths" & Res[,2] %in% genome,3] <- paste(GenomeBuild,"_ChrLengths.txt",sep="")
        if(length(grep(".gtf:",GenePos)) > 0){
          GenePosToUpdate <- Res[Res[,1] %in% "Gene Positions" & Res[,2] %in% genome,3] <- paste(gtf,"_PCK.bed",sep="")
        }
     }
     
      AllSections <- unique(Res[,1])
      file.rename(ConfigOut,gsub(".ini","preupdatedanno.ini",ConfigOut))
      file.create(ConfigOut,showWarnings=T)
      for(i in 1:length(AllSections)){
        write.table(paste("[",AllSections[i],"]",sep=""),ConfigOut,col.names=F,row.names=F,quote=F,append=T)
        VariablesInSection <- as.vector(Res[Res[,1] %in% AllSections[i],2])
        for(j in 1:length(VariablesInSection)){
             Variable <- VariablesInSection[j]
             Value  <- as.vector(Res[Res[,1] %in% AllSections[i] & Res[,2] %in% VariablesInSection[j],3])
             write.table(paste(Variable,Value,sep=" = "),ConfigOut,col.names=F,row.names=F,quote=F,append=T)      
        }
        write.table("",ConfigOut,col.names=F,row.names=F,quote=F,append=T)      
      }

     
     
}

