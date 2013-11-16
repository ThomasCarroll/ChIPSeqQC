library(DiffBind)
library(DiffBindCRI)
system("mkdir -p DiffBind")
sampleSheetPreDBA <- read.delim("SampleSheet.csv",sep=",",h=T)
sampleSheetPreDBA <- sampleSheetPreDBA[is.na(sampleSheetPreDBA[,"toMerge"]),]

WkgDir <- getwd()

SampleID <- sampleSheetPreDBA[,"SampleName"]
Tissue <- sampleSheetPreDBA[,"Tissue"]
Factor <- sampleSheetPreDBA[,"Factor"]
Condition <- sampleSheetPreDBA[,"Condition_1"]
Replicate <- sampleSheetPreDBA[,"Replicate"]
bamReads <- file.path("bamFiles",sampleSheetPreDBA[,"Processed_bamFileName"])

bamControl <- vector("character",nrow(sampleSheetPreDBA))
ControlID <- vector("character",nrow(sampleSheetPreDBA))
for(i in 1:length(bamControl)){
if(any(sampleSheetPreDBA[,"SampleName"] %in% sampleSheetPreDBA[i,"InputToUse"])){
  bamControl[i] <-  file.path("bamFiles",as.vector(sampleSheetPreDBA[sampleSheetPreDBA[,"SampleName"] %in% sampleSheetPreDBA[i,"InputToUse"],"Processed_bamFileName"]))
  ControlID[i] <- as.vector(sampleSheetPreDBA[sampleSheetPreDBA[,"SampleName"] %in% sampleSheetPreDBA[i,"InputToUse"],"SampleName"])
}
else if(any(sampleSheetPreDBA[,"GenomicsID"] %in% sampleSheetPreDBA[i,"InputToUse"])){
  bamControl[i] <-  file.path("bamFiles",as.vector(sampleSheetPreDBA[sampleSheetPreDBA[,"GenomicsID"] %in% sampleSheetPreDBA[i,"InputToUse"],"Processed_bamFileName"]))
    ControlID[i] <- as.vector(sampleSheetPreDBA[sampleSheetPreDBA[,"GenomicsID"] %in% sampleSheetPreDBA[i,"InputToUse"],"SampleName"])
}else{
   bamControl[i] <- NA
   ControlID[i] <- NA
}
}

Peaks <- sampleSheetPreDBA[,"Macs_name"]

#Peaks <- rep("/lustre/mib-cri/carrol09/Work/Kirsten_From_Rory/Peak.bed",length=nrow(sampleSheetPreDBA))

print(sessionInfo())

PeakCaller <- rep("macs",nrow(sampleSheetPreDBA))
PeakCaller[is.na(Peaks)] <- NA
sampleSheet <- data.frame(SampleID,Tissue,Factor,Condition,Replicate,bamReads,bamControl,ControlID,Peaks)
sampleSheet <- sampleSheet[!is.na(sampleSheet[,"Peaks"]),]
write.table(sampleSheet,"forDBAsampleSheet.csv",sep=",",quote=F,row.names=F)
DBAObj <- dba(sampleSheet="forDBAsampleSheet.csv",caller="macs")

png(file.path(WkgDir,"DiffBind","Occupancy_Heatmap.png"))
plot(DBAObj)
dev.off()

DBAObj$config$parallelPackage <- 2
DBAObj <-  dba.count(DBAObj)

png(file.path(WkgDir,"DiffBind","Affinity_Heatmap.png"))
plot(DBAObj)
dev.off()

png(file.path(WkgDir,"DiffBind","Affinity_PCA.png"))
dba.plotPCA(DBAObj)
dev.off()
