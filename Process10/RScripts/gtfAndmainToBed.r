Args <- commandArgs(trailingOnly = TRUE)

BedOut <- Args[3]
#Args[1] <- "/bio11array1/carrol09/Test5/TestRun/annotation/mm9_Ensembl/Mus_musculus.NCBIM37.67.gtf.bed"
#Args[2] <- "/bio11array1/carrol09/Test5/TestRun/annotation/mm9_Ensembl/mmusculus_gene_ensembl__transcript__main.txt"
Temp <- read.delim(Args[1],sep="\t",h=F)
Temp2 <- read.delim(Args[2],sep="\t",h=F)[,c(62,64,116,122,164)]
#Temp <- read.delim("/bio11array1/carrol09/Test4/Output.bed",sep="\t",h=F)
#Temp2 <- read.delim("/bio11array1/carrol09/Test4/annotation/mm10_Ensembl/mmusculus_gene_ensembl__transcript__main.txt",sep="\t",h=F)[,c(61,63,115,121,163)]
Temp3 <- Temp2[Temp2[,2] %in% "protein_coding" &  Temp2[,4] %in% "KNOWN",]
Temp4 <- merge(Temp,Temp3,by.x=4,by.y=1,all.x=F,all.y=T)
AllGeneNames <-  unique(Temp4[,16])
Temp5 <- NULL
for(i in 1:length(AllGeneNames)){
    AllPCKTranscriptsForGene <- Temp4[Temp4[,16] %in% AllGeneNames[i],]
    AllPCKTranscriptsForGene <- AllPCKTranscriptsForGene[!is.na(as.vector(AllPCKTranscriptsForGene[,2])) & !is.na(AllPCKTranscriptsForGene[,3]) & !is.na(AllPCKTranscriptsForGene[,4]),]
    CombinedTranscript <-   cbind(unique(as.vector(AllPCKTranscriptsForGene[,2]))[1],min(AllPCKTranscriptsForGene[,3]),max(AllPCKTranscriptsForGene[,4]),unique(as.vector(AllPCKTranscriptsForGene[,6]))[1],as.vector(AllGeneNames[i]))
    Temp5 <- rbind(Temp5,CombinedTranscript)
}
colnames(Temp5) <- c("Chr","Start","End","Strand","Symbol")
write.table(Temp5,BedOut,row.names=F,quote=F,sep="\t")
