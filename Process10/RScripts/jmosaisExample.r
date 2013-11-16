bsub -q bioinformatics -J R_8G -R "rusage[mem=16000]" /lustre/mib-cri/carrol09/QCPaper/R-3.0.1/bin/Rscript mosaics.r /lustre/mib-cri/carrol09/QCPaper/bamFiles/IL-LNCAP-LM-EV-BICALUTAMIDE-INPUT-D7.bwa_Processed.bed 150
bsub -q bioinformatics -J R_8G -R "rusage[mem=16000]" /lustre/mib-cri/carrol09/QCPaper/R-3.0.1/bin/Rscript mosaics.r /lustre/mib-cri/carrol09/QCPaper/bamFiles/S37-AR-CHIP-LNCAP-LM-EV-BICALUTAMIDE.bwa_Processed.bed 150
bsub -q bioinformatics -J R_8G -R "rusage[mem=16000]" /lustre/mib-cri/carrol09/QCPaper/R-3.0.1/bin/Rscript mosaics.r /lustre/mib-cri/carrol09/QCPaper/bamFiles/S38-AR-CHIP-LNCAP-LM-EV-BICALUTAMIDE.bwa_Processed.bed 150

library(jmosaics)

bin1_noToy <- readBins(type = c("chip","input"),fileName = c("/lustre/mib-cri/carrol09/QCPaper/S37-AR-CHIP-LNCAP-LM-EV-BICALUTAMIDE.bwa_Processed.bed_fragL150_bin200.txt","/lustre/mib-cri/carrol09/QCPaper/IL-LNCAP-LM-EV-BICALUTAMIDE-INPUT-D7.bwa_Processed.bed_fragL150_bin200.txt"))
bin2_noToy <- readBins(type = c("chip","input"),fileName = c("/lustre/mib-cri/carrol09/QCPaper/S38-AR-CHIP-LNCAP-LM-EV-BICALUTAMIDE.bwa_Processed.bed_fragL150_bin200.txt","/lustre/mib-cri/carrol09/QCPaper/IL-LNCAP-LM-EV-BICALUTAMIDE-INPUT-D7.bwa_Processed.bed_fragL150_bin200.txt"))

origin_bin_toy=list(bin1_noToy,bin2_noToy)
bin<- readBinsMultiple(origin_bin)
fit1 <- mosaicsFit(bin[[1]], analysisType = "IO", bgEst="automatic")
fit2 <- mosaicsFit(bin[[2]], analysisType = "IO", bgEst="automatic")
fit <- list(fit1, fit2)
result<-jmosaicsPattern(fit, region_length=1, FDR=0.05, thres=c(10,10),type=c('B', 'E', 'Pattern'), patternInfo='TRUE')
