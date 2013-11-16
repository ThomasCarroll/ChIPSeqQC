optionsList <- vector("list",length=8)
optionsList["PeakDirectory"] <- "/lustre/mib-cri/carrol09/Work/acrapplace/DB Sites/"
optionsList["RankOrder"] <- "Rank"
optionsList["ColumnIndexForRank"] <- 6
optionsList["config"] <-"/lustre/mib-cri/stark01/20130610_MohammedH_JC_PRChIPs/Xenografts/Config/config.ini"

#optionsList["PeakFileNameColumn"] <- "/lustre/mib-cri/stark01/20130610_MohammedH_JC_PRChIPs/Xenografts/"

CallMotifsCheck <- function(Caller,WkgDir=getwd(),Config="Config"){
   if(Caller != "Custom"){
   Flags <- GetFlags(WkgDir,Config)
   LogicalCall <- slot(Flags,paste("Call",Caller,"Motifs",sep="")) == "Yes"
   return(LogicalCall)
   }else{
       return(TRUE)
   }
}

WkgDir <- getwd()

GetPipelinebase <- function(WkgDir=getwd(),ConfigDirectory="Config"){
  require(raster)
  ConfigToRead = file.path(WkgDir,ConfigDirectory,"config.ini")
  ConfigFile <- readIniFile(ConfigToRead)
  
  PipelineBase <- ConfigFile[ConfigFile[,2] %in% "BaseLocation",3]
  return(PipelineBase)
}  


RunPeakMotifingPipeline <- function(SampleSheet,WkgDir=WkgDir,JobString,MaxJobs=75,PLs=PipelineLocations,Config="Config",Caller="Macs",NPeaks=500,optionsList=Null){
PipelineBase <- GetPipelinebase(dirname(unlist(optionsList["config"])),"")

source(file.path(PipelineBase,"/RScripts/Workflow_Functions3.r"))


#if(CallMotifsCheck(Caller,WkgDir,Config)){
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
  if(Caller != "Custom" || is.null(optionsList)){
#      optionsList["Name"] <- "/lustre/mib-cri/carrol09/Work/acrapplace/"
      



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

  }
  
#    Pipeline = "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/MotifProcessPipeline2.xml"
    PipeName <- paste("SS_",Caller,"Motifs",sep="")
    SSPeakProfile_WfMeta <- MakeWorkFlowXML2(JobString,Variables,Specialisations,Pipeline,PipeName,WkgDir,Config,maxJobs=75)
    saveXML(SSPeakProfile_WfMeta,file=file.path(PLs@WorkFlowDir,paste("SS_",Caller,"MotifProcess.xml",sep="")))
    system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,paste("SS_",Caller,"MotifProcess.xml",sep="")),sep=""),wait=T)

  }else if(Caller == "Custom" & !is.null(optionsList)){

  WkgDir <- getwd()
   print("Here")
  PeakDirectory <- unlist(optionsList[unlist("PeakDirectory")])
  PeakFileNameColumn <- optionsList["ColumnIndexForName"]
  ColumnIndexForRank <- optionsList["ColumnIndexForRank"]
  RankOrder <- optionsList["RankOrder"]
  #PeaksToRunDir <- dir(optionsList["Peaks"])
  Peaks <- dir(PeakDirectory)
  
  Pipeline <- getPipelinesPath("motifpipeline",dirname(unlist(optionsList["config"])),"")
  workflowExec <- getWorkflowParam("Executable",dirname(unlist(optionsList["config"])),"")
  javaExec <- getExecPath("java",dirname(unlist(optionsList["config"])),"")
  mode <- tolower(getWorkflowParam("Mode",dirname(unlist(optionsList["config"])),""))
  PipelineBase <- GetPipelinebase(dirname(unlist(optionsList["config"])),"")
  pythonExec <- getExecPath("python",dirname(unlist(optionsList["config"])),"")
  pipelineRun <- paste(javaExec," -jar ",workflowExec," --mode=",mode,sep="")
  rexec <- getExecPath("rexec",dirname(unlist(optionsList["config"])),"")
  meme <- getExecPath("meme",dirname(unlist(optionsList["config"])),"")
  ame <- getExecPath("ame",dirname(unlist(optionsList["config"])),"")
  rlibs <- getLibraryPath("rlibs",dirname(unlist(optionsList["config"])),"")
  perllibs <- getLibraryPath("perllibs",dirname(unlist(optionsList["config"])),"")
  pythonlibs <- getLibraryPath("pythonlibs",dirname(unlist(optionsList["config"])),"")

  genomeChrLengths <- GetChrLengthsFromConfig(dirname(unlist(optionsList["config"])),"")
  print("Here3")
  PLs <- GetImportantLocations(dirname(unlist(optionsList["config"])),"")

  BamDir <- PLs@BamDir
  MemeFormatdb <- GetTFDB(dirname(unlist(optionsList["config"])),"")
  genome <- GetGenomeFromConfig(dirname(unlist(optionsList["config"])),"")
  genomeChrLengths <- GetChrLengthsFromConfig(dirname(unlist(optionsList["config"])),"")
  fastaFromGenome <- GetGenomeBuildFromConfig(dirname(unlist(optionsList["config"])),"")
  Variables  <- c(getwd(),file.path(PeakDirectory,""),file.path(BamDir,""),genome,genomeChrLengths,fastaFromGenome,NPeaks,RankOrder,ColumnIndexForRank,MemeFormatdb,ame,meme,rexec,PipelineBase,pythonExec,rlibs,perllibs,pythonlibs,genomeChrLengths)
  names(Variables) <- c("WorkingDirectory","PeakDirectory","BamDirectory","genome","genomeFile","Fasta","NPeaks","RankOrder","RankColumn","MotifDatabaseLocation","ame","meme","rexec","pipelineBase","python","rlibs","perllibs","pythonlibs","chrlengths")

  print("Here2")
  Specialisations  <- vector("list",length=length(Peaks))
  for(i in 1:length(Specialisations)){
    Specialisations[[i]] <- c("",Peaks[i],gsub("\\.bed","",Peaks[i]))
      names(Specialisations[[i]]) <- c("SamplePeakDir","PeakFileName","Test")
  }
  names(Specialisations) <-  gsub("\\.bed","",Peaks)

  #    Pipeline = "/lustre/mib-cri/carrol09/Work/MyPipe/Process10/src/main/pipelines/MotifProcessPipeline2.xml"
  PipeName <- paste("SS_",Caller,"Motifs",sep="")
  SSPeakProfile_WfMeta <- MakeWorkFlowXML2("Custom",Variables,Specialisations,Pipeline,PipeName,dirname(dirname(unlist(optionsList["config"]))),"Config",maxJobs=75)
  saveXML(SSPeakProfile_WfMeta,file=file.path(PLs@WorkFlowDir,paste("SS_",Caller,"MotifProcess.xml",sep="")))
  system(paste(pipelineRun," ",file.path(PLs@WorkFlowDir,paste("SS_",Caller,"MotifProcess.xml",sep="")),sep=""),wait=T)


  }
  

#}

}



optionsList <- list()
optionsList["PeakDirectory"] <- "/lustre/mib-cri/stark01/20130610_MohammedH_JC_PRChIPs/SitesForMotifAnalysis/"
optionsList["RankOrder"] <- "Rank"
optionsList["ColumnIndexForRank"] <- 5
optionsList["ColumnIndexForName"] <- 4
optionsList["config"] <-"/lustre/mib-cri/stark01/20130610_MohammedH_JC_PRChIPs/Xenografts/Config/config.ini"

RunPeakMotifingPipeline(SampleSheet="",WkgDir=getwd(),JobString="Custom",MaxJobs=75,PLs="",Config="Config",Caller="Custom",NPeaks=500,optionsList=optionsList)


