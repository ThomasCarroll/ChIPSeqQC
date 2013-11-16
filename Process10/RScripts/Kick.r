#!/lustre/mib-cri/carrol09/Work/MyPipe/R/R-2.15.0/bin/Rscript --vanilla
RLIBSVar = "/bio11array1/carrol09/Test4/lib/R/library/" 

suppressPackageStartupMessages(library("optparse",lib=RLIBSVar))
#library("methods")
suppressPackageStartupMessages(library("sp",lib=RLIBSVar))
suppressPackageStartupMessages(library("raster",lib=RLIBSVar))
suppressPackageStartupMessages(library("XML",lib=RLIBSVar))
suppressPackageStartupMessages(library("RJSONIO",lib=RLIBSVar))

#print("All Libraries loaded")
## Create the directory structure based on config
system(paste("mkdir -p Config",sep=""),wait=T)

##Get the base!
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
PipeBase <- paste(as.character(unlist(strsplit(script.name,.Platform$file.sep))[-c(length(unlist(strsplit(script.name,.Platform$file.sep)))-1,length(unlist(strsplit(script.name,.Platform$file.sep))))]),collapse=.Platform$file.sep)
print(paste("Running scripts from ",PipeBase,sep=""))
## Source Workflow functions 
source(file.path(PipeBase,"RScripts","Workflow_Functions3.r"))
#
WkgDir <- getwd()
#CreateDirStruct(WkgDir,"Config",PipeLineLocations)
#print(script.name)
 
ConfigOut <- file.path(WkgDir,"Config","config.ini")


option_list <- list(


  make_option(c("-c","--config"),type="character",help="The pipeline config file",default=file.path(PipeBase,"Config/config.ini")),
	make_option(c("-g","--genome"),type="character",help="The desired genome to be used"),

	make_option(c("--workingDirectory"),type="character",help="The desired genome to be used",default=""),
	make_option(c("--bamDirectory"),type="character",help="The desired genome to be used",default=""),
	make_option(c("--fastqDirectory"),type="character",help="The desired genome to be used",default=""),

	make_option(c("--addSLXIDs"),type="character",help="The desired genome to be used",default=""),
	make_option(c("--addSampleNames"),type="character",help="The desired genome to be used",default=""),
	make_option(c("--addProjects"),type="character",help="The desired genome to be used",default=""),

	make_option(c("--callMacsAll"),type="character",help="Call,QC and find motifs in Macs Peaks",default=""),
	make_option(c("--callMacsPeaks"),type="character",help="Call Macs Peaks",default=""),
	make_option(c("--callMacsMotifs"),type="character",help="Call Macs Motif",default=""),
	make_option(c("--callMacsPeakProfile"),type="character",help="Perform QC on Macs Peaks",default=""),
	make_option(c("--callMacsBetweenPeaks"),type="character",help="Perform QC between Macs Peaks",default=""),
  	
	make_option(c("--callSicerAll"),type="character",help="Call,QC and find motifs in Sicer Peaks",default=""),	
	make_option(c("--callSicerPeaks"),type="character",help="Call Sicer Peaks",default=""),
	make_option(c("--callSicerMotifs"),type="character",help="Call Sicer Motif",default=""),
	make_option(c("--callSicerPeakProfile"),type="character",help="Perform QC on Sicer Peaks",default=""),
	make_option(c("--callSicerBetweenPeaks"),type="character",help="Perform QC between Sicer Peaks",default=""),
	

	make_option(c("--callTPICsAll"),type="character",help="Call,QC and find motifs in TPICs Peaks",default=""),	
	make_option(c("--callTPICsPeaks"),type="character",help="Call TPICs Peaks",default=""),
	make_option(c("--callTPICsMotifs"),type="character",help="Call TPICs Motif",default=""),
	make_option(c("--callTPICsPeakProfile"),type="character",help="Perform QC on TPICs Peaks",default=""),
	make_option(c("--callTPICsBetweenPeaks"),type="character",help="Perform QC between TPICs Peaks",default=""),

	make_option(c("--AutoMerging"),type="character",help="Merge samples by samplename",default=""),	
	make_option(c("--JustDummy"),type="character",help="Makes dummy samplesheet,directory structure and config",default="")	

)

opt <- parse_args(OptionParser(option_list=option_list,description="\tWelcome to the CRUK Cambridge Institute ChIP-seq Pipeline\n\tWritten by Tom Carroll\n\tSuraj Menon and Rory Stark\n\tthomas.carroll@cruk.cam.ac.uk",prog="ChIPseqPipeline"))
if(!file.exists(ConfigOut)){
  Res <- readIniFile(opt$config)
##Match Options and Config
}else{
  Res <- readIniFile(ConfigOut)
}

namesVersusConfig <- tolower(names(opt))
for(i in 1:length(opt)){
  if(namesVersusConfig[i] %in% Res[,2]){
    if(opt[[i]] != ""){
      Res[Res[,2] %in% namesVersusConfig[i],3] <- opt[[i]]
    }
  }
}

Res[Res[,2] %in% "workingdirectory" & Res[,3] %in% "",3] <- file.path(WkgDir,"")
Res[Res[,2] %in% "bamdirectory" & Res[,3] %in% "",3] <- file.path(WkgDir,"bamFiles","")
Res[Res[,2] %in% "fastqdirectory" & Res[,3] %in% "",3] <- file.path(WkgDir,"FQFiles","")
Res[Res[,2] %in% "tempdirectory" & Res[,3] %in% "",3] <- file.path(WkgDir,"Temp","")
Res[Res[,2] %in% "locationsdirectory" & Res[,3] %in% "",3] <- file.path(WkgDir,"Locations","")
Res[Res[,2] %in% "workflowdirectory" & Res[,3] %in% "",3] <- file.path(WkgDir,"Workflow","")
Res[Res[,2] %in% "fraglengthdirectory" & Res[,3] %in% "",3] <- file.path(WkgDir,"Fragment_Lengths","")
Res[Res[,2] %in% "coveragedirectory" & Res[,3] %in% "",3] <- file.path(WkgDir,"Coverage","")

Res[Res[,2] %in% "macsdirectory" & Res[,3] %in% "",3] <- file.path(WkgDir,"Peaks","Macs","")
Res[Res[,2] %in% "sicerdirectory" & Res[,3] %in% "",3] <- file.path(WkgDir,"Peaks","Sicer","")
Res[Res[,2] %in% "tpicsdirectory" & Res[,3] %in% "",3] <- file.path(WkgDir,"Peaks","TPICs","")


#print(Res)
#if(Res[Res[,2] %in% "genome",3] == ""){
#  stop("Need to specify a Genome!")
  #exit()
#}


Res[Res[,2] %in% "taskdirectories",3] <- file.path(PipeBase,"src","main",basename(Res[Res[,2] %in% "taskdirectories",3]))

AllSections <- unique(Res[,1])[!unique(Res[,1]) %in% c("Custom Scripts","PipeLine_Base","Pipelines")]
file.rename(ConfigOut,gsub(".ini","old.ini",ConfigOut))
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

ScriptSections <- unique(Res[,1])[unique(Res[,1]) %in% "Custom Scripts"]
for(i in 1:length(ScriptSections)){
  write.table(paste("[",ScriptSections[i],"]",sep=""),ConfigOut,col.names=F,row.names=F,quote=F,append=T)
  VariablesInSection <- as.vector(Res[Res[,1] %in% ScriptSections[i],2])
  for(j in 1:length(VariablesInSection)){
       Variable <- VariablesInSection[j]
       ValueWithoutBase  <- as.vector(Res[Res[,1] %in% ScriptSections[i] & Res[,2] %in% VariablesInSection[j],3])
       Temp <- unlist(strsplit(ValueWithoutBase,.Platform$file.sep))
       Value <- file.path(PipeBase,paste(Temp[(grep("Process10",Temp)[1]+1):length(Temp)],collapse=.Platform$file.sep))
       write.table(paste(Variable,Value,sep=" = "),ConfigOut,col.names=F,row.names=F,quote=F,append=T)      
  }
  write.table("",ConfigOut,col.names=F,row.names=F,quote=F,append=T)      
}


pipeSections <- unique(Res[,1])[unique(Res[,1]) %in% "Pipelines"]
for(i in 1:length(pipeSections)){
  write.table(paste("[",pipeSections[i],"]",sep=""),ConfigOut,col.names=F,row.names=F,quote=F,append=T)
  VariablesInSection <- as.vector(Res[Res[,1] %in% pipeSections[i],2])
  for(j in 1:length(VariablesInSection)){
       Variable <- VariablesInSection[j]
       ValueWithoutBase  <- basename(as.vector(Res[Res[,1] %in% pipeSections[i] & Res[,2] %in% VariablesInSection[j],3]))
       Value <- file.path(PipeBase,"src","main","pipelines",ValueWithoutBase)
      # print(Value)
       write.table(paste(Variable,Value,sep=" = "),ConfigOut,col.names=F,row.names=F,quote=F,append=T)      
  }
  write.table("",ConfigOut,col.names=F,row.names=F,quote=F,append=T)      
}

write.table(paste("[","PipeLine_Base","]",sep=""),ConfigOut,col.names=F,row.names=F,quote=F,append=T)
write.table(paste("BaseLocation",PipeBase,sep=" = "),ConfigOut,col.names=F,row.names=F,quote=F,append=T)


system(paste("mkdir -p ",Res[Res[,2] %in% "locationsdirectory",3],sep=""),wait=T) 
#Res[Res[,2] %in% "locationsdirectory",3]

SLXIDsFile <- file.path(Res[Res[,2] %in% "locationsdirectory",3],"SLXIDsFromLims.txt")
ProjectsFile <- file.path(Res[Res[,2] %in% "locationsdirectory",3],"ProjectsFromLims.txt")
SampleNamesFile <- file.path(Res[Res[,2] %in% "locationsdirectory",3],"SamplenamesFromLims.txt")

#print(SLXIDsFile)

addSLXIDs <- opt$addSLXIDs
if(addSLXIDs != ""){
  SLXIDs <- as.vector(unlist(strsplit(addSLXIDs,",")))
  if(file.exists(SLXIDsFile) & file.info(SLXIDsFile)$size > 0){
#  print("mate")
    TempSLXIDs <- as.vector(read.delim(SLXIDsFile,sep="",h=F)[,1])
    write.table(unique(c(TempSLXIDs,SLXIDs)),SLXIDsFile,col.names=F,row.names=F,quote=F)
  }else{
    write.table(unique(c(SLXIDs)),SLXIDsFile,col.names=F,row.names=F,quote=F)
  }
}

addProjects <- opt$addProjects
if(addProjects != ""){
  Projects <- as.vector(unlist(strsplit(addProjects,",")))
}else{
  Projects <- ""
}  

if(file.exists(ProjectsFile) & file.info(ProjectsFile)$size > 0){
  TempProjects <- as.vector(read.delim(ProjectsFile,sep="",h=F)[,1])
  write.table(unique(c(basename(getwd()),TempProjects,Projects)),ProjectsFile,col.names=F,row.names=F,quote=F)
}else{
  write.table(unique(c(basename(getwd()),Projects)),ProjectsFile,col.names=F,row.names=F,quote=F)
}


addSampleNames <- opt$addSampleNames
if(addSampleNames != ""){
  SampleNames <- as.vector(unlist(strsplit(addSampleNames,",")))
  if(file.exists(SampleNamesFile) & file.info(SampleNamesFile)$size > 0){
    TempSampleNames <- as.vector(read.delim(SampleNamesFile,sep="",h=F)[,1])
    write.table(unique(c(TempSampleNames,SampleNames)),SampleNamesFile,col.names=F,row.names=F,quote=F)
  }else{
    write.table(unique(c(SampleNames)),SampleNamesFile,col.names=F,row.names=F,quote=F)
  }
}
  if(!file.exists(SLXIDsFile)){
      file.create(SLXIDsFile,showWarnings = F)
  }
  if(!file.exists(ProjectsFile)){
      file.create(ProjectsFile,showWarnings = F)
  }
    if(!file.exists(SampleNamesFile)){
      file.create(SampleNamesFile,showWarnings = F)
  }
  
#WkgDir <- Args[2]
getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
JobString <- getRandString()
## Parse from config important locations

PipeLineLocations <- GetImportantLocations(WkgDir,"Config")


## Create the directory structure based on config
CreateDirStruct(WkgDir,"Config",PipeLineLocations)

## Make sure samplesheet is unlocked before starting.
UnlockSampleSheet(WkgDir)
if(opt$JustDummy != "Yes"){
cat("Launching main pipeline\n.\n.\n.\n.\n.\n.\n.\n.\n.\n.\n.\n.\n.\n")
RunMainPipeline(WkgDir,JobString,MaxJobs=75,PipeLineLocations,"Config")
}else{
SampleSheet <- GetSampleSheetOrDummy(WkgDir,"SampleSheet.csv")
write.table(SampleSheet,file.path(WkgDir,"SampleSheet.csv"),sep=",",quote=F,row.names=F)
}
write.table("Complete",file.path(PipeLineLocations@WorkFlowDir,paste(JobString,"_Main_0.txt",sep="")),col.names=T,row.names=F,sep=",",quote=F)

