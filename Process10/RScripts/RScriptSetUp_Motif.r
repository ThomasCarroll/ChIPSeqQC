getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
source("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/RScripts/Workflow_Functions.r")

Args <- commandArgs(trailingOnly = TRUE)
library(raster)
library(XML)
#ConfigFile <- readIniFile(file.path(Args[2],"config.ini"))
ConfigFile <- readIniFile(file.path(getwd(),"Temp","config.ini"))
LocationsDir <- ConfigFile[ConfigFile[,2] %in% "tempdirectory",3]
WkgDir<- ConfigFile[ConfigFile[,2] %in% "workingdirectory",3]
BamDir<- ConfigFile[ConfigFile[,2] %in% "bamdirectory",3]
Genome <- ConfigFile[ConfigFile[,2] %in% "genome",3]
GenomeFileOptions <-  ConfigFile[ConfigFile[,1] %in% "Genomes",]
GenomeFile <- GenomeFileOptions[GenomeFileOptions[,2] %in% tolower(Genome),3]



Config <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/Config/Config.txt",sep="\t",header=F)


CallMotifs <- ConfigFile[ConfigFile[,2] %in% "callmacsmotifs",3]

if(CallMotifs %in% "Yes"){
	CallMotifs <- TRUE
}else{
	CallMotifs <- FALSE	
}

if(CallMotifs){

sampleSheet <- read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=",")
MacsPeaks <- sampleSheet[!is.na(sampleSheet[,"Macs_name"]),"Macs_name"]
SampleNames <- sampleSheet[!is.na(sampleSheet[,"Macs_name"]),"GenomicsID"]
MacsPeaks <- gsub(".xls",".bed",MacsPeaks)






TempFull <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/MotifMeta.xml",header=F,stringsAsFactors=F,quote="")
TempFull[grep("RunDirectory_ToChange",TempFull[,1]),1] <- gsub("RunDirectory_ToChange",getwd(),TempFull[grep("RunDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("BamDirectory_ToChange",TempFull[,1]),1] <- gsub("BamDirectory_ToChange",BamDir,TempFull[grep("BamDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1] <- gsub("WorkingDirectory_ToChange",WkgDir,TempFull[grep("WorkingDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("Genome_ToChange",TempFull[,1]),1] <- gsub("Genome_ToChange",GenomeFile,TempFull[grep("Genome_ToChange",TempFull[,1]),1])
TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1] <- gsub("TempDirectory_ToChange",LocationsDir,TempFull[grep("TempDirectory_ToChange",TempFull[,1]),1])
TempFull[grep("MotifDatabaseLocation_ToChange",TempFull[,1]),1] <- gsub("MotifDatabaseLocation_ToChange","/lustre/mib-cri/carrol09/Work/MyPipe/NewMeme/TranfacMatrix.meme",TempFull[grep("MotifDatabaseLocation_ToChange",TempFull[,1]),1])


SpecialisationSet <- TempFull[grep("specialisation identifier",TempFull[,1]):max(grep("specialisation>",TempFull[,1])),1]
AllSets <- rep(list(SpecialisationSet),length(MacsPeaks))
BigSet <- vector("character")
for(i in 1:length(AllSets)){
  AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])] <- gsub("Indentifier_ToChange",getRandString(len=6),AllSets[[i]][grep("Indentifier_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("Test_ToChange",AllSets[[i]])] <- gsub("Test_ToChange",SampleNames[i],AllSets[[i]][grep("Test_ToChange",AllSets[[i]])])
  AllSets[[i]][grep("PeakFileLocation_ToChange",AllSets[[i]])] <- gsub("PeakFileLocation_ToChange",MacsPeaks[i],AllSets[[i]][grep("PeakFileLocation_ToChange",AllSets[[i]])])
  BigSet <- c(BigSet,AllSets[[i]])
}

FirstAllExceptSpecialisations <- TempFull[c(1:grep("specialisations",TempFull[,1])[1]),1]
SecondAllExceptSpecialisations <- TempFull[c(grep("specialisations",TempFull[,1])[2]:nrow(TempFull)),1]
Tommy <- c(FirstAllExceptSpecialisations,BigSet,SecondAllExceptSpecialisations)
write.table(Tommy,file=file.path(LocationsDir,"RunMotif.xml"),qmethod="escape",quote=F,sep="\t",col.names=F,row.names=F)
cat("Submitting jobs!............")
system(paste("java -Xmx1G -jar /lustre/mib-cri/carrol09/MyPipe/workflow-all-1.2-SNAPSHOT.jar --mode=lsf ",file.path(LocationsDir,"RunMotif.xml"),sep=""),wait=TRUE,intern=FALSE)



#sampleSheet <- read.delim(file.path(WkgDir,"SampleSheet.csv"),sep=",",stringsAsFactors=F)
sampleSheet <- ReadAndLock(file.path(WkgDir,"SampleSheet.csv"),WkdDir,SAF=F,napTime=5)

for(i in 1:nrow(sampleSheet)){
	SampleToLookFor <- sampleSheet[i,"GenomicsID"]
	KnownMotifFile <- dir(path=file.path(WkgDir,"Motif",SampleToLookFor,"Known"),pattern="*me.txt$",full.names=T)
	if(all(c(length(KnownMotifFile) > 0,file.info(KnownMotifFile)$size > 0))){
		sampleSheet[i,"Known_Motif_File"] <- KnownMotifFile
		print(KnownMotifFile)
		try(
		DataIn <- read.delim(KnownMotifFile,sep=" ",comment.char="#",skip=11,header=F)
		,silent=T)
		try(
		TempSignificantTable <- cbind(sampleSheet[i,"GenomicsID"],as.vector(DataIn[,6]),as.vector(DataIn[,7]),as.vector(DataIn[,11]),gsub(")","",as.vector(DataIn[,14])))
		,silent=T)
		try(
		sampleSheet[i,"Known_Significant_Motifs"]  <- nrow(DataIn)
		,silent=T)
		try(
		if(exists(as.character(bquote(SignificantTable)))){
			SignificantTable <- rbind(SignificantTable,TempSignificantTable)
		}else{
			SignificantTable <- TempSignificantTable
		}
		,silent=T)
		
	}
}
colnames(SignificantTable) <- c("Genomics_ID","Motif_ID","Motif_Name","Motif_Fisher_PValue","Motif_Fisher_AdjPValue")
write.table(SignificantTable,file.path(WkgDir,"Motif","Summary_Significance_Table_KnownMotifs.txt"),sep="\t")

MatchToMatrix <- read.delim("/lustre/mib-cri/carrol09/Work/MyPipe/Motifs/matrix_list.txt",sep="\t",header=F)

for(i in 1:nrow(sampleSheet)){
	SampleToLookFor <- as.vector(sampleSheet[i,"GenomicsID"])
	print(SampleToLookFor)
	DremeMotifFile <- dir(path=file.path(WkgDir,"Motif",SampleToLookFor,"Denovo","dreme_out"),pattern="*.xml$",full.names=T)
	if(all(c(length(DremeMotifFile) > 0,file.info(DremeMotifFile)$size > 0))){
		doc = xmlTreeParse(DremeMotifFile, useInternal = TRUE)
		top = xmlRoot(doc)
		NOfMotifs <- length(names(top[["motifs"]]))
		for(k in 1:NOfMotifs){
#			ProbLine <- top[["motifs"]][[k]][["match"]]
			LineAtrr <- xmlAttrs(top[["motifs"]][[k]])
			TempSignificantTableDreme <- matrix(c(SampleToLookFor,LineAtrr[1:9]),nrow=1,byrow=T)
			if(exists(as.character(bquote(SignificantTableDreme)))){
				SignificantTableDreme <- rbind(SignificantTableDreme,TempSignificantTableDreme)
			}else{
				SignificantTableDreme <- TempSignificantTableDreme
			}
		}
		sampleSheet[i,"Dreme_Motif_File"] <- dir(path=file.path(WkgDir,"Motif",SampleToLookFor,"Denovo","dreme_out"),pattern="*.html$",full.names=T)
		sampleSheet[i,"Dreme_Significant_Motifs"] <- NOfMotifs
	}
}



for(i in 1:nrow(sampleSheet)){
	SampleToLookFor <- as.vector(sampleSheet[i,"GenomicsID"])
	DremeTomTomFile <- dir(path=file.path(WkgDir,"Motif",SampleToLookFor,"Denovo","dreme_tomtom_out"),pattern="*.txt$",full.names=T)
	if(all(c(length(DremeTomTomFile) > 0,file.info(DremeTomTomFile)$size > 0))){
		TempTomTom <- read.delim(DremeTomTomFile,sep="\t",header=T)
		if(nrow(TempTomTom) > 0){
		PerSampleTomTom <- cbind(SampleToLookFor,TempTomTom)
		if(exists(as.character(bquote(SignificantTomTomDreme)))){
			SignificantTomTomDreme <- rbind(SignificantTomTomDreme,PerSampleTomTom)
		}else{
			SignificantTomTomDreme <- PerSampleTomTom 
		}
			sampleSheet[i,"Dreme_Significant_Motifs_Matched_To_Known"] <- length(unique(PerSampleTomTom[,2]))
		}else{
			sampleSheet[i,"Dreme_Significant_Motifs_Matched_To_Known"] <- 0
		}
	}
}

PerSampleMotifsDremeTomTom <- paste(SignificantTomTomDreme[,1],SignificantTomTomDreme[,2],sep="__")
PerSampleMotifsDremeOnly <- paste(SignificantTableDreme[,1],SignificantTableDreme[,3],sep="__")

for(i in 1:length(PerSampleMotifsDremeTomTom)){
	Temp <- matrix(SignificantTableDreme[PerSampleMotifsDremeOnly %in% PerSampleMotifsDremeTomTom[i],],nrow=1,byrow=T)
	TempBiggerDremeTomTom <- cbind(SignificantTomTomDreme[i,],Temp)
	if(exists(as.character(bquote(BiggerDremeTomTom)))){
		BiggerDremeTomTom <- rbind(BiggerDremeTomTom,TempBiggerDremeTomTom)
	}else{
		BiggerDremeTomTom <- TempBiggerDremeTomTom
	}
}
	
NewDremeTomTom <- BiggerDremeTomTom[,-c(12,13,14,15)]
#colnames(NewDremeTomTom)[12:17] <-  c("nsites","p","n","P-Value","E-Value","Unerased_Evalue")
DremeTomTom <- merge(NewDremeTomTom,MatchToMatrix,by.x=3,by.y=1,all.x=T,all.y=F)
DremeTomTom <- cbind(DremeTomTom[,c(2)],DremeTomTom[,c(1,18,19)],DremeTomTom[,-c(1,2,18,19)])
colnames(DremeTomTom)[c(1,3,4,14,15,16,17,18,19)] <- c("SampleName","Tranfac_ID","Transcription_Factor","nsites","p","n","P-Value","E-Value","Unerased_Evalue")
write.table(DremeTomTom,file=file.path(WkgDir,"Motif","Summary_Significance_Table_Dreme_Denovo_Motifs.txt"),row.names=F,sep=",")

for(i in 1:nrow(sampleSheet)){
	SampleToLookFor <- sampleSheet[i,"GenomicsID"]
	MemeMotifFile <- dir(path=file.path(WkgDir,"Motif",SampleToLookFor,"Denovo","meme_out"),pattern="*.xml$",full.names=T)
	if(all(c(length(MemeMotifFile) > 0,file.info(MemeMotifFile)$size > 0))){
		doc = xmlTreeParse(MemeMotifFile, useInternal = TRUE)
		top = xmlRoot(doc)
		NOfMotifs <- length(names(top[["motifs"]]))
		for(k in 1:NOfMotifs){
			LineAtrr <- xmlAttrs(top[["motifs"]][[k]])
			TempSignificantTablememe <- matrix(c(SampleToLookFor,LineAtrr[1:9]),nrow=1,byrow=T)
			if(exists(as.character(bquote(SignificantTablememe)))){
				SignificantTablememe <- rbind(SignificantTablememe,TempSignificantTablememe)
			}else{
				SignificantTablememe <- TempSignificantTablememe
			}
		}
	
		sampleSheet[i,"meme_Motif_File"] <- dir(path=file.path(WkgDir,"Motif",SampleToLookFor,"Denovo","meme_out"),pattern="*.html$",full.names=T)
		sampleSheet[i,"meme_Significant_Motifs"] <- NOfMotifs	
	}
}


for(i in 1:nrow(sampleSheet)){
	SampleToLookFor <- sampleSheet[i,"GenomicsID"]
	memeTomTomFile <- dir(path=file.path(WkgDir,"Motif",SampleToLookFor,"Denovo","meme_tomtom_out"),pattern="*.txt$",full.names=T)
	if(all(c(length(memeTomTomFile) > 0,file.info(memeTomTomFile)$size > 0))){
		TempTomTom <- read.delim(memeTomTomFile,sep="\t",header=T)
		if(nrow(TempTomTom) > 0){
		PerSampleTomTom <- cbind(SampleToLookFor,TempTomTom)
		if(exists(as.character(bquote(SignificantTomTommeme)))){
			SignificantTomTommeme <- rbind(SignificantTomTommeme,PerSampleTomTom)
		}else{
			SignificantTomTommeme <- PerSampleTomTom 
		}
		sampleSheet[i,"meme_Significant_Motifs_Matched_To_Known"] <- length(unique(PerSampleTomTom[,2]))
		}else{
		sampleSheet[i,"meme_Significant_Motifs_Matched_To_Known"] <- 0
		}
	}
}

PerSampleMotifsmemeTomTom <- paste(SignificantTomTommeme[,1],SignificantTomTommeme[,2],sep="__")
PerSampleMotifsmemeOnly <- paste(SignificantTablememe[,1],SignificantTablememe[,3],sep="__")

for(i in 1:length(PerSampleMotifsmemeTomTom)){
	Temp <- matrix(SignificantTablememe[PerSampleMotifsmemeOnly %in% PerSampleMotifsmemeTomTom[i],],nrow=1,byrow=T)
	TempBiggermemeTomTom <- cbind(SignificantTomTommeme[i,],Temp)
	if(exists(as.character(bquote(BiggermemeTomTom)))){
		BiggermemeTomTom <- rbind(BiggermemeTomTom,TempBiggermemeTomTom)
	}else{
		BiggermemeTomTom <- TempBiggermemeTomTom
	}
}

#             id,name,width,sites,ic,re,llr,e_value,bayes_threshold,elapsed_time

NewmemeTomTom <- BiggermemeTomTom[,-c(12,13,14)]
#colnames(NewmemeTomTom)[12:17] <-  c("nsites","p","n","P-Value","E-Value","Unerased_Evalue")
memeTomTom <- merge(NewmemeTomTom,MatchToMatrix,by.x=3,by.y=1,all.x=T,all.y=F)
memeTomTom <- cbind(memeTomTom[,c(2)],memeTomTom[,c(1,19,20)],memeTomTom[,-c(1,2,19,20)])
colnames(memeTomTom)[c(1,3,4,14,15,16,17,18,19,20)] <- c("SampleName","Tranfac_ID","Transcription_Factor","width","sites","ic","re","llr","e_value","bayes_threshold")
write.table(memeTomTom,file=file.path(WkgDir,"Motif","Summary_Significance_Table_meme_Denovo_Motifs.txt"),row.names=F,sep=",")


#write.table(sampleSheet,file=file.path(WkgDir,"SampleSheet.csv"),row.names=F,sep=",")
WriteAndUnlock(sampleSheet,file.path(WkgDir,"SampleSheet.csv"))	

}else{cat("No Motifs To Call")}
write.table("Complete",file.path(LocationsDir,paste(Args[1],"_MainMotifProcess.txt",sep="")),col.names=T,row.names=F,sep=",",quote=F)
