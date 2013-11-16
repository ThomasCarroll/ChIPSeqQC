library(DiffBind)
library(DiffBindCRI)
library(raster)
system("mkdir -p DiffBind")

print(sessionInfo())

Args <- commandArgs(trailingOnly = TRUE)


WkgDir <- getwd()
sampleSheetPreDBA <- read.delim("SampleSheet.csv",sep=",",h=T)
sampleSheetPreDBA <- sampleSheetPreDBA[is.na(sampleSheetPreDBA[,"toMerge"]) & sampleSheetPreDBA[,"Analysis_State"] %in% "RunMe",]
colnames(sampleSheetPreDBA)
ConfigFile <- readIniFile(file.path(getwd(),"Config","config.ini"))

caller <- Args[1]

diffbind.config = function(config,param,default) {
   retval  = NULL
   if(sum(config[,1] %in% "DiffBind")>0) {
     DiffBindOptions = config[config[,1] %in% "DiffBind",]
      retval <- DiffBindOptions[DiffBindOptions[,2] %in% param,3]
      if(length(retval)==0) {
         retval=NULL   	
      }   	
   }
   if(is.null(retval)) {
      retval=default
   }
   return(retval)	
}

maxpeaks   = diffbind.config(ConfigFile,"max_peaks", 100000)
minoverlap = diffbind.config(ConfigFile,"min_overlap", 0.2)
minmembers = diffbind.config(ConfigFile,"min_members", 3)
threshold  = as.numeric(diffbind.config(ConfigFile,"threshold", 0.1))
doedgeR    = diffbind.config(ConfigFile,"edgeR", TRUE)
dodeseq    = diffbind.config(ConfigFile,"DESeq", FALSE)

methodlist = NULL
if(doedgeR) methodlist = c(methodlist,DBA_EDGER)
if(dodeseq) methodlist = c(methodlist,DBA_DESEQ)

peaks = dba.pipeline(bCorPlot=F)

diffbind.analysis = function(peaks,minoverlap,maxpeaks,methodlist,threshold,caller,bCount=T,bAnalyze=T) {

   dba.save(peaks,sprintf('peaks_%s',caller),file.path(WkgDir,"DiffBind"))

   png(file.path(WkgDir,"DiffBind",sprintf("Occupancy_CorHeatmap_%s.png",caller)),width=1024,height=1024)
   plot(peaks,attributes=c(DBA_ID,DBA_CALLER))
   dev.off()

   rates = dba.overlap(peaks,mode=DBA_OLAP_RATE)
   olap = which(rates <= maxpeaks)[1]
   olap = max(minoverlap,olap)
   png(file.path(WkgDir,"DiffBind",sprintf("Overlap_Rate_%s.png",caller)),width=1024,height=800)
   plot(rates,type="b",xlab="# Samples",ylab="# Peaks")
   points(olap,rates[olap],type='p',pch=19,col="red")
   dev.off()

   if(bCount) {
      counts <-  dba.count(peaks,minOverlap=olap,bCorPlot=F)
      dba.save(counts,sprintf('counts_%s',caller),file.path(WkgDir,"DiffBind"))
   } else {
      counts = dba.load(sprintf('counts_%s',caller),file.path(WkgDir,"DiffBind"))	
   }

   png(file.path(WkgDir,"DiffBind",sprintf("Affinity_CorHeatmap_%s.png",caller)),width=1024,height=1024)
   plot(counts)
   dev.off()

   png(file.path(WkgDir,"DiffBind",sprintf("Affinity_PCA_%s.png",caller)),width=800,height=800)
   dba.plotPCA(counts,dotSize=2)
   dev.off()

   counts    = dba.contrast(counts,minMembers=minmembers)
   contrasts = dba.show(counts,bContrasts=T)

   if(!is.null(contrasts)) {
   	
   	  if(bAnalyze) {
         counts = dba.analyze(counts,method=methodlist,bCorPlot=F)
         dba.save(counts,sprintf('analysis_%s',caller),file.path(WkgDir,"DiffBind"))
      } else {
      	 counts = dba.load(sprintf('analysis_%s',caller),file.path(WkgDir,"DiffBind"))
      }
      
      contrasts = dba.show(counts,bContrasts=T)
      write.csv(contrasts,row.names=F, file=file.path(WkgDir,"DiffBind",sprintf("Contrasts_%s.csv",caller)))	
   
      for(i in 1:nrow(contrasts)) {
         for(method in methodlist) {
      	    
      	    if(method==DBA_EDGER) methodname = "edgeR"
      	    if(method==DBA_DESEQ) methodname = "DESeq"
      	    
            png(file.path(WkgDir,"DiffBind",sprintf("Contrast%d_MA_%s_%s.png",i,methodname,caller)),width=1200,height=1024)
            dba.plotMA(counts,method=method,contrast=1)
            dev.off()

            rep1 = dba.report(counts,method=method,contrast=i,th=1,DataType=DBA_DATA_FRAME)
            write.csv(rep1,row.names=F, file=file.path(WkgDir,"DiffBind",sprintf("Contrast%d_Full_Report_%s_%s.csv",i,methodname,caller)))
            
            if(sum(rep1$FDR<threshold)>0) { 
                
              png(file.path(WkgDir,"DiffBind",sprintf("Contrast%d_CorHeatmap_%s_%s.png",i,methodname,caller)),width=1024,height=1024)
              dba.plotHeatmap(counts,method=method,contrast=i)
              dev.off()
              
              png(file.path(WkgDir,"DiffBind",sprintf("Contrast%d_Heatmap_%s_%s.png",i,methodname,caller)),width=1024,height=1024)
              dba.plotHeatmap(counts,method=method,contrast=i,correlations=F)
              dev.off()
              
              png(file.path(WkgDir,"DiffBind",sprintf("Contrast%d_PCA_%s_%s.png",i,methodname,caller)),width=800,height=800)
              dba.plotPCA(counts,method=method,contrast=i,dotSize=2)
              dev.off()

              rep2 = dba.report(counts,method=method,contrast=i,th=threshold,DataType=DBA_DATA_FRAME)
              write.csv(rep2,row.names=F, file=file.path(WkgDir,"DiffBind",sprintf("Contrast%d_DB_Report_%s_%s.csv",i,methodname,caller)))
            }           
            #diffbind_contrastReport(counts,i,methodname,caller,rep1,rep2,rate,olap)     
         }
      }
   }
   
   #diffbind_report()
   
   return(list(DBA=counts, rates=rates,rate=olap,contrasts=contrasts))
}

#peaks$config$parallelPackage <- 1

#analysis       = diffbind.analysis(peaks,minoverlap,maxpeaks,methodlist,threshold,"all")
if(tolower(caller) == "macs"){
peaks.macs     = dba(peaks,peaks$masks$macs)
analysis.macs  = diffbind.analysis(peaks.macs,minoverlap,maxpeaks,methodlist,threshold,"macs")
}
if(tolower(caller) == "sicer"){
peaks.sicer    = dba(peaks,peaks$masks$sicer)
analysis.sicer = diffbind.analysis(peaks.sicer,minoverlap,maxpeaks,methodlist,threshold,"sicer")
}
if(tolower(caller) == "tpics"){
peaks.tpic     = dba(peaks,peaks$masks$tpic)
analysis.tpic  = diffbind.analysis(peaks.tpic,minoverlap,maxpeaks,methodlist,threshold,"tpic")
}
