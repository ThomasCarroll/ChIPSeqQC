getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
WkgDir <- getwd()

Args <- commandArgs(trailingOnly = TRUE)
library(raster)
UserName = system("id -nu",wait=TRUE,intern=TRUE)
MegaName <- paste("cri.camres.org\\\\",UserName,"@uk-cri-larc01",sep="")

ConfigFile <- readIniFile(file.path(Args[2],"config.ini"))
ConfigFile <- readIniFile(file.path(WkgDir,"Temp","config.ini"))

#LocationsDir <- file.path(WkgDir,"config.ini")

LocationsDir <- ConfigFile[ConfigFile[,2] %in% "tempdirectory",3]
WkgDir<- ConfigFile[ConfigFile[,2] %in% "workingdirectory",3]
if(file.info(file.path(LocationsDir,"Lims_Info.txt"))$size != 0 & file.info(file.path(LocationsDir,"ActualLocations.txt"))$size != 0){ 
MajorInfo <- read.delim(file.path(LocationsDir,"Lims_Info.txt"),sep="\t",header=F)
RemainingInfo <- matrix(nrow=nrow(MajorInfo),ncol=22,data=NA)
SampleSheet <- cbind(MajorInfo,RemainingInfo)


MinorInfo <- gsub(".*=> ","",as.vector(read.delim(file.path(LocationsDir,"ActualLocations.txt"),sep="\t",header=F)[,1]))
MinoMat <- matrix(MinorInfo,ncol=7,byrow=T)
MinoMatFinal <- cbind(matrix(unlist(strsplit(MinoMat[,2],".bwa")),ncol=2,byrow=T)[,1],paste(MinoMat[,7],":",MinoMat[,4],"/",MinoMat[,2],sep=""))


SampleSheet[,5] <- MinoMatFinal[match(SampleSheet[,1],MinoMatFinal[,1]),2]
SampleSheet[is.na(SampleSheet[,5]),5] <- "No_Lims_Location_Known"

SampleSheet[,5] <- gsub("archive.crnet.org",MegaName,as.vector(SampleSheet[,5]),fixed=T)

#colnames(SampleSheet) <- c("GenomicsID","SampleName","Run","Lane","Location","Tissue","Factor","Antibody","Condition_1","Condition_2","Replicate","InputToUse","toMerge","bamFileName","Processed_bamFileName","Original","delRand","Excluded","Filtered","Unique","DuplicationRate","InputFileNameToUse","derivedFrom","Macs_Fragment_Length","Sissr_Fragment_Length","Correlation_Fragment_Length","Coverage_Fragment_Length","MacsPeaks","Macs_name","SicerPeaks","Sicer_Name","TPICsPeaks","TPICS_Name")
colnames(SampleSheet) <- c("GenomicsID","SampleName","Run","Lane","Location","Tissue","Factor","Antibody","Condition_1","Condition_2","Replicate","InputToUse","toMerge","bamFileName","Processed_bamFileName","Original","delRand","Excluded","Filtered","Unique","DuplicationRate","InputFileNameToUse","derivedFrom","Macs_Fragment_Length","Sissr_Fragment_Length","Correlation_Fragment_Length","Coverage_Fragment_Length")

if(file.exists("SampleSheet.csv")){

Oldssht <- read.table(file.path(WkgDir,"SampleSheet.csv"),sep=",",header=T)
if(ncol(Oldssht) - ncol(SampleSheet) != 0 ){
ExtraForOld <- ncol(Oldssht) - ncol(SampleSheet)
ExpandOld <- cbind(SampleSheet,matrix(nrow=nrow(SampleSheet),ncol=ExtraForOld))
colnames(ExpandOld)[ncol(SampleSheet):ncol(ExpandOld)] <- colnames(Oldssht)[ncol(SampleSheet):ncol(ExpandOld)]
NewSheet <- rbind(Oldssht,ExpandOld)
}else{
NewSheet <- rbind(Oldssht,SampleSheet)
}

SampleSheet <- NewSheet[match(unique(NewSheet[,1]),NewSheet[,1]),]
}


LocationToGet <-  gsub("/solexa/","/",as.vector(SampleSheet[!as.vector(SampleSheet[,"Location"]) %in% "No_Lims_Location_Known" & !is.na(as.vector(SampleSheet[,"Location"])),"Location"]))
BamNames <-  as.vector(SampleSheet[!as.vector(SampleSheet[,"Location"]) %in% "No_Lims_Location_Known"  & !is.na(as.vector(SampleSheet[,"Location"])),"GenomicsID"])

#if(length(LocationToGet) > 0){

TempFull <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/GetBamsMeta.xml",header=F,stringsAsFactors=F,quote="")
TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1] <- gsub("WorkingDirectory_ToChange",WkgDir,TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1] <- gsub("TempDirectory_ToChange",LocationsDir,TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1])


SpecialisationSet <- TempFull[grep("specialisation identifier",TempFull[,1]):max(grep("specialisation>",TempFull[,1])),1]
AllSets <- rep(list(SpecialisationSet),length(LocationToGet))
BigSet <- vector("character")
for(i in 1:length(AllSets)){
  AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])] <- gsub("Indentifier_ToChange",getRandString(len=6),AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("Bam_To_Get",AllSets[[i]])] <- gsub("Bam_To_Get",LocationToGet[i],AllSets[[i]][grep("Bam_To_Get",AllSets[[i]])],fixed=T)
  AllSets[[i]][grep("Bam_Name_To_Get",AllSets[[i]])] <- gsub("Bam_Name_To_Get",BamNames[i],AllSets[[i]][grep("Bam_Name_To_Get",AllSets[[i]])])
  BigSet <- c(BigSet,AllSets[[i]])
}

FirstAllExceptSpecialisations <- TempFull[c(1:grep("specialisations",TempFull[,1])[1]),1]
SecondAllExceptSpecialisations <- TempFull[c(grep("specialisations",TempFull[,1])[2]:nrow(TempFull)),1]
Tommy <- c(FirstAllExceptSpecialisations,BigSet,SecondAllExceptSpecialisations)
write.table(Tommy,file=file.path(LocationsDir,"TrialBamGet.xml"),qmethod="escape",quote=F,sep="\t",col.names=F,row.names=F)
cat("Submitting jobs!............")
system(paste("java -jar /lustre/mib-cri/carrol09/MyPipe/workflow-all-1.2-SNAPSHOT.jar --mode=lsf ",file.path(LocationsDir,"TrialBamGet.xml"),sep=""),wait=TRUE,intern=FALSE)
cat("Jobs Submitted!\n")


BamFilesReceived <- dir(path=file.path(WkgDir,"bamFiles"),pattern="*.bam$")
if(any(grep("Processed",BamFilesReceived))){BamFilesReceived <- BamFilesReceived[-grep("Processed",BamFilesReceived)]}
if(any(grep("Realign",BamFilesReceived))){BamFilesReceived <- BamFilesReceived[-grep("Realign",BamFilesReceived)]}


SampleSheet <- as.matrix(SampleSheet)
BamID <- matrix(as.vector(unlist(strsplit(BamFilesReceived,"\\.bwa"))),ncol=2,byrow=T)
 
for(i in 1:nrow(SampleSheet)){
	if(!any(c(grep("Processed",SampleSheet[i,"bamFileName"]),grep("Realign",SampleSheet[i,"bamFileName"])))){
		if(any(BamID %in% SampleSheet[i,"GenomicsID"])){
			SampleSheet[i,"bamFileName"] <- BamFilesReceived[BamID %in% SampleSheet[i,"GenomicsID"]]
		}
	}
}

if(file.exists(file.path(LocationsDir,"metadata.txt"))){
	MetaSampleSheets <- read.delim(file.path(LocationsDir,"metadata.txt"),stringsAsFactors=F,header=F)[,1]
	for(i in 1:length(MetaSampleSheets)){
		MetaInfo <- read.delim(MetaSampleSheets[i],stringsAsFactors=F,header=T,sep=",")
		for(k in 1:nrow(SampleSheet)){
			GenomicID <- SampleSheet[k,"GenomicsID"]
			if(any(MetaInfo[,1] %in% GenomicID)){
				if(any(colnames(MetaInfo) %in% "Condition")){
					MetaInfoFull <- MetaInfo[MetaInfo[,1] %in% GenomicID,c("Tissue","Factor","Antibody","Condition","Replicate","InputToUse","toMerge")]
					SampleSheet[k,c("Tissue","Factor","Antibody","Condition_1","Replicate","InputToUse","toMerge")] <- as.vector(as.matrix(MetaInfoFull))
				}else{
					MetaInfoFull <- MetaInfo[MetaInfo[,1] %in% GenomicID,c("Tissue","Factor","Antibody","Condition_1","Condition_2","Replicate","InputToUse","toMerge")]
					SampleSheet[k,c("Tissue","Factor","Antibody","Condition_1","Condition_2","Replicate","InputToUse","toMerge")] <- as.vector(as.matrix(MetaInfoFull))
				}
				
			}
		}
		
	}
	
}
write.table(SampleSheet,file.path(WkgDir,"SampleSheet.csv"),col.names=T,row.names=F,sep=",",quote=F)
}
write.table("Complete",file.path(LocationsDir,paste(Args[1],"_MainBamGet.txt",sep="")),col.names=T,row.names=F,sep=",",quote=F)
