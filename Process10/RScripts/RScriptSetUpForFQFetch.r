getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))

Args <- commandArgs(trailingOnly = TRUE)
library(raster)
ConfigFile <- readIniFile(file.path(Args[2],"config.ini"))
ConfigFile <- readIniFile(file.path(getwd(),"Temp","config.ini"))
LocationsDir <- ConfigFile[ConfigFile[,2] %in% "tempdirectory",3]
WkgDir<- ConfigFile[ConfigFile[,2] %in% "workingdirectory",3]
BamDir<- ConfigFile[ConfigFile[,2] %in% "bamdirectory",3]
 
UserName = system("id -nu",wait=TRUE,intern=TRUE)
MegaName <- paste("cri.camres.org\\\\",UserName,"@uk-cri-larc01",sep="")



if(file.exists(file.path(LocationsDir,"ActualFQLocations.txt")) & file.info(file.path(LocationsDir,"ActualFQLocations.txt"))$size > 0){

	Stuff <- as.vector(read.delim(file.path(LocationsDir,"ActualFQLocations.txt"),sep="\t",header=F)[,1])
	IDIndicies <- grep("[[:space:]]ID",Stuff)
	IDIndicies <- c(IDIndicies,length(Stuff))
	StartIndex <- IDIndicies[1:(length(IDIndicies)-1)]
	EndIndex <- IDIndicies[2:(length(IDIndicies))]

	MinoMat <- matrix(nrow=1,ncol=8)	

	for(i in 1:length(StartIndex)){
		if(i < length(StartIndex)){
			TempStuff <- Stuff[(StartIndex[i]+1):(EndIndex[i]-1)]
		}else{
			TempStuff <- Stuff[(StartIndex[i]+1):(EndIndex[i])]
		}
		TempSLXID <- Stuff[(StartIndex[i])]
		TempMinorInfo <- gsub(".*=> ","",TempStuff)
		if(length(TempMinorInfo) > 7){
			print(i)
		}
		TempMinoMat <- matrix(TempMinorInfo,ncol=7,byrow=T)
		TempMinoMat <- cbind(rep(gsub(".*=> ","",TempSLXID),nrow(TempMinoMat)),TempMinoMat)
		for(j in 1:nrow(TempMinoMat)){
			RunID <- strsplit(strsplit(TempMinoMat[j,][5],"Runs")[[1]][2],"_")[[1]][3]
			SeqID <- strsplit(TempMinoMat[j,3],"_se")[[1]][1] 
			TempMinoMat[j,2] <- paste(gsub(".*=> ","",TempSLXID),RunID,SeqID,sep=".")
		}
		MinoMat <- rbind(MinoMat,TempMinoMat)
	}
	MinoMat <- matrix(MinoMat,byrow=F,ncol=8)
	MinoMat <- MinoMat[-1,]
	MinoMatFinal <- cbind(MinoMat[,2],paste(MinoMat[,8],":",MinoMat[,5],"/",MinoMat[,3],sep=""))




MinoMatFinal[,2] <- gsub("archive.crnet.org",MegaName,as.vector(MinoMatFinal[,2]),fixed=T)

write.table(MinoMatFinal,file.path(LocationsDir,"FastQPositions.txt"),sep="\t",row.names=F,quote=F)


TempFull <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/GetFQMeta.xml",header=F,stringsAsFactors=F,quote="")
TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1] <- gsub("WorkingDirectory_ToChange",getwd(),TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1] <- gsub("TempDirectory_ToChange",LocationsDir,TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1])

SpecialisationSet <- TempFull[grep("specialisation identifier",TempFull[,1]):max(grep("specialisation>",TempFull[,1])),1]
AllSets <- rep(list(SpecialisationSet),nrow(MinoMatFinal))

BigSet <- vector("character")
for(i in 1:length(AllSets)){
  AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])] <- gsub("Indentifier_ToChange",getRandString(len=6),AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("Bam_To_Get",AllSets[[i]])] <- gsub("Bam_To_Get",MinoMatFinal[i,2],AllSets[[i]][grep("Bam_To_Get",AllSets[[i]])],fixed=T)
  AllSets[[i]][grep("Bam_Name_To_Get",AllSets[[i]])] <- gsub("Bam_Name_To_Get",MinoMatFinal[i,1],AllSets[[i]][grep("Bam_Name_To_Get",AllSets[[i]])])
  BigSet <- c(BigSet,AllSets[[i]])
}

FirstAllExceptSpecialisations <- TempFull[c(1:grep("specialisations",TempFull[,1])[1]),1]
SecondAllExceptSpecialisations <- TempFull[c(grep("specialisations",TempFull[,1])[2]:nrow(TempFull)),1]
Tommy <- c(FirstAllExceptSpecialisations,BigSet,SecondAllExceptSpecialisations)
write.table(Tommy,file=file.path(LocationsDir,"TrialFQGet.xml"),qmethod="escape",quote=F,sep="\t",col.names=F,row.names=F)
cat("Submitting jobs!............")
system(paste("java -jar /lustre/mib-cri/carrol09/MyPipe/workflow-all-1.2-SNAPSHOT.jar --mode=lsf ",file.path(LocationsDir,"TrialFQGet.xml"),sep=""),wait=TRUE,intern=FALSE)
cat("Jobs Submitted!\n")
}
write.table("Complete",file.path(LocationsDir,paste(Args[1],"_MainFQProcess.txt",sep="")),col.names=T,row.names=F,sep=",",quote=F)


