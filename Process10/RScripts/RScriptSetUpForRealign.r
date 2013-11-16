getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))

Args <- commandArgs(trailingOnly = TRUE)
library(raster)
ConfigFile <- readIniFile(file.path(Args[2],"config.ini"))
#ConfigFile <- readIniFile(file.path(getwd(),"Temp","config.ini"))

LocationsDir <- ConfigFile[ConfigFile[,2] %in% "tempdirectory",3]
WkgDir<- ConfigFile[ConfigFile[,2] %in% "workingdirectory",3]
BamDir<- ConfigFile[ConfigFile[,2] %in% "bamdirectory",3]
  
#Config <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/Config/Config.txt",sep="\t",header=F)
GENOME_PARAMETER <- ConfigFile[ConfigFile[,2] %in% "genome",3]


TempFull <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/MultiGenomeGetter.xml",header=F,stringsAsFactors=F,quote="")
TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1] <- gsub("WorkingDirectory_ToChange",getwd(),TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1] <- gsub("TempDirectory_ToChange",LocationsDir,TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("BamDirectory_ToChange",TempFull[,1]),1] <- gsub("BamDirectory_ToChange",BamDir,TempFull[grep("BamDirectory_ToChange",TempFull[,1]),1])


TempFull[grep("genome_ToChange",TempFull[,1]),1] <- gsub("genome_ToChange",getwd(),TempFull[grep("genome_ToChange",TempFull[,1]),1])

sampleSheet <- read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=",")
Samples <- gsub("\\.bam","",sampleSheet[!is.na(sampleSheet[,"bamFileName"]) & ! sampleSheet[,"Location"] %in% "No_Lims_Location_Known","bamFileName"])


SpecialisationSet <- TempFull[grep("specialisation identifier",TempFull[,1]):max(grep("specialisation>",TempFull[,1])),1]
AllSets <- rep(list(SpecialisationSet),length(Samples))

BigSet <- vector("character")
for(i in 1:length(AllSets)){
  AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])] <- gsub("Indentifier_ToChange",getRandString(len=6),AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("BamName_ToChange",AllSets[[i]])] <- gsub("BamName_ToChange",Samples[i],AllSets[[i]][grep("BamName_ToChange",AllSets[[i]])])
  BigSet <- c(BigSet,AllSets[[i]])
}



FirstAllExceptSpecialisations <- TempFull[c(1:grep("specialisations",TempFull[,1])[1]),1]
SecondAllExceptSpecialisations <- TempFull[c(grep("specialisations",TempFull[,1])[2]:nrow(TempFull)),1]
Tommy <- c(FirstAllExceptSpecialisations,BigSet,SecondAllExceptSpecialisations)
write.table(Tommy,file=file.path(LocationsDir,"GetGenome.xml"),qmethod="escape",quote=F,sep="\t",col.names=F,row.names=F)



cat("Submitting jobs!............")
system(paste("java -jar /lustre/mib-cri/carrol09/MyPipe/workflow-all-1.2-SNAPSHOT.jar --mode=lsf ",file.path(LocationsDir,"GetGenome.xml"),sep=""),wait=TRUE,intern=FALSE)
cat("Jobs Submitted!\n")
files <- dir(path=BamDir,pattern="*.info",full.names = T)
print(files)
Genomes <- vector("character",length=length(files))
for(i in 1:length(files)){
	Genomes[i] <- read.delim(files[i],header=F,sep="\t")
}


Genomes <- as.vector(unlist(Genomes))
fileNames <- gsub(".info","",basename(files))
GenomeMat <- cbind(paste(fileNames,".bam",sep=""),Genomes)
write.table(cbind(paste(fileNames,".bam",sep=""),Genomes),file.path(WkgDir,"bamFiles","GenomeVersions.txt"),sep="\t",row.names=F,col.names=F,quote=F)


GenomesToAlign <- GenomeMat[!tolower(GenomeMat[,2]) %in% tolower(GENOME_PARAMETER),1]
if(file.exists(file.path(LocationsDir,"SLXToAlign.txt"))){
	unlink(file.path(LocationsDir,"SLXToAlign.txt"))
}
if(length(GenomesToAlign) > 0){
SLXToAlign <- vector("character",length=length(GenomesToAlign))
for(i in 1:length(GenomesToAlign)){
if(any(sampleSheet[,"bamFileName"] %in% GenomesToAlign[i])){
	SLXToAlign[i] <- strsplit(as.vector(sampleSheet[sampleSheet[,"bamFileName"] %in% GenomesToAlign[i],1]),"\\.")[[1]][1]
}	
}
SLXToAlign <- SLXToAlign[SLXToAlign != ""]
write.table(SLXToAlign,file.path(LocationsDir,"SLXToAlign.txt"),sep="\t",row.names=F,col.names=F,quote=F)
}
write.table("Complete",file.path(LocationsDir,paste(Args[1],"_MainGenomeProcess.txt",sep="")),col.names=T,row.names=F,sep=",",quote=F)

