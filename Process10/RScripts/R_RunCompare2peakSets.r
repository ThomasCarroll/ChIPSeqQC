Args <- commandArgs(trailingOnly = TRUE)
library(GenomicRanges)

WkgDir <- getwd()

source("/lustre/mib-cri/carrol09/Work/MyPipe/Process10/RScripts/Analysis_Functions.r")

#Args[1] <- "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/Trial/TestForTest/Peaks/Sicer/SLX-5967.MDAMB134POLIIPHOPHOS2B.985.s_7.bwa.homo_sapiens_Processed/SLX-5967.MDAMB134POLIIPHOPHOS2B.985.s_7.bwa.homo_sapiens_Processed-W200-G600-islands-summary-FDR0.01"
#Args[2] <- "/lustre/mib-cri/carrol09/Work/PipelinePracticeSet/Trial/TestForTest/Peaks/Sicer/SLX-5607.908.s_7.bwa.homo_sapiens_Processed/SLX-5607.908.s_7.bwa.homo_sapiens_Processed-W200-G600-islands-summary-FDR0.01"
#Args[3] <- "FirstOutFile.txt"
#Args[4] <- "ResultsFile.txt"

print(Args[1])
print(Args[2])
print(Args[3])
print(Args[4])


GRanges_1 <-  GetGRanges(Args[1],simple=T)
GRanges_2 <-  GetGRanges(Args[2],simple=T)



GRList <- GRangesList(GRanges_1,GRanges_2)
names(GRList) <- c(basename(Args[1]),basename(Args[2]))

GR_Reps <- CountRepPeaks(GRList)

#print("GotHere")
#GR_Reps

#print(GRanges_1)
#print(GRanges_2)

JS <- JaccardScore(GRanges_1,GRanges_2)
#JS <- as.vector(JaccardOverGRangeslist(GRList)[1,3])
#print("GotHere2")
AllOverLaps <- elementMetadata(GR_Reps)$PeakInRep_Num
AllPeaks <- length(AllOverLaps)
InAtLeast1 <- length(AllOverLaps[AllOverLaps > 1])
JSofPeaks <-  InAtLeast1/AllPeaks
#print("GotHere3")
LengthOfSample1 <- length(GR_Reps[elementMetadata(GR_Reps)[,1] > 0])
LengthOfSample2 <- length(GR_Reps[elementMetadata(GR_Reps)[,2] > 0])

PercentOf1in2 <-  InAtLeast1/LengthOfSample1
PercentOf2in1 <-  InAtLeast1/LengthOfSample2

AllRes <- c(LengthOfSample1,LengthOfSample2,AllPeaks,InAtLeast1,PercentOf1in2,PercentOf2in1,JS,JSofPeaks)
ScoreNames <- c("Peaks_Sample1","Peaks_Sample2","Union_Peaks_Of_Samples","Intersection_Peaks_Of_Samples","Percent_Of_Intersection_Peaks_In_Sample1","Percent_Of_Intersection_Peaks_In_Sample2","Jaccard_Score_Of_Peaks_Overlap","Jaccard_Score_Of_BasePair_Overlap")
ResMatrix <- cbind(ScoreNames,AllRes)
write.table(ResMatrix,Args[4],quote=F,row.names=F,col.names=F,sep="\t")
write.table(as.data.frame(GR_Reps),Args[3],quote=F,row.names=F,col.names=T,sep="\t")


