Args <- commandArgs(trailingOnly = TRUE)
library(mosaics)
bedFile <- Args[1]
fragLen <- Args[2]
Bins <- constructBins(infile=bedFile,fileFormat="bed", outfileLoc="/bamFiles/",byChr=FALSE, useChrfile=FALSE, chrfile=NULL, excludeChr=NULL,PET=FALSE, fragLen=fragLen, binSize=200, capping=0 )
