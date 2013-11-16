library(rGADEM)
library(Biostrings)
Arguments <- commandArgs(trailingOnly = T)
getRandString<-function(len=12) return(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE),collapse=''))
Fasta <- Arguments[1]
FastaSet <- read.DNAStringSet(Fasta, format="fasta",use.names=TRUE)
gadem <-GADEM(FastaSet,verbose=1,genome=Hsapiens)
save(gadem,file=paste(".fa","_GademResults.RData",Fasta))