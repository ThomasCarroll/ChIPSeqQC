ReadAndLock <- function(ss,WkdDir,SAF=T,napTime=5){
  if(file.exists(gsub(".csv",".LOCK",ss))){
    while(file.exists(gsub(".csv",".LOCK",ss))){
      Sys.sleep(napTime)
    }
    write.table("Locked",gsub(".csv",".LOCK",ss))
    sampleSheet <- read.delim(ss,stringsAsFactors=SAF,sep=",")
  }else{
    write.table("Lcoked",gsub(".csv",".LOCK",ss))
    sampleSheet <- read.delim(ss,stringsAsFactors=SAF,sep=",")
  }
  return(sampleSheet)

}

WriteAndUnlock <- function(sampleSheet,ss){
   if(file.exists(gsub(".csv",".LOCK",ss))){
     write.table(sampleSheet,ss,sep=",",row.names=F,quote=F)
     unlink(gsub(".csv",".LOCK",ss))
   }else{
     write.table(sampleSheet,ss,sep=",",row.names=F,quote=F)
   }
}

