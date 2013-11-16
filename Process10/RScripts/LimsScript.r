findSLXIDS <- function(SLXID=NULL,LimsOutPut="LimsInfo.txt"){
library(XML)
library(RCurl)
library(RJSONIO)
#SLXID <- c("SLX-7483")

LimsLocations <- NULL
if(!is.null(SLXID)){
  for(i in 1:length(SLXID)){
    new <- getURL(paste("http://genomicsequencing.cruk.cam.ac.uk:8080/glsintapi/runsContainingLibraries?slxId=",SLXID[i],sep=""))
    Runs <- xmlToList(new)
    slxidRuns <- vector("character")
    slxidTypes <- vector("character")

    TempCount <- which(Runs[1,] != "MiSeq Run")

    for(k in 1:length(TempCount)){
       slxidRuns[k] <- Runs[,TempCount[k]]$runFolder
       slxidTypes[k] <- Runs[,TempCount[k]]$runType
    }
    Temp <- cbind(rep(SLXID[i],length(slxidRuns)),slxidRuns,slxidTypes)
    #Temp <- Temp[-1,]
    colnames(Temp)  <- c("SLXID","Runs","Types")

     TempFull <- Temp

    for(j in 1:nrow(TempFull)){
       Files <- xmlToList(getURL(paste("http://genomicsequencing.cruk.cam.ac.uk:8080/glsintapi/fullDetailsOfRun?runId=",TempFull[j,"Runs"],sep="")))
       Lanes <- Files[rownames(Files) == "flowcell",][[1]][-c(1,length(Files[rownames(Files) == "flowcell",][[1]]))]
       incount = 0
       Temp2Full <- NULL
       for(l in 1:length(Lanes)){
         if(Lanes[[l]]$slxId == SLXID[i]){
            incount = incount+1
            print(l)
            JustSamples <- Lanes[[l]][names(Lanes[[l]]) == "sample"]
            vecSampleLocations <-  vector("character")
            vecSampleNames <-  vector("character")
            if(any(names(unlist(JustSamples)) == "sample.file.url")){
              for(f in 1:length(JustSamples)){
                JustSamples[[f]]$name
                vecSample <- matrix(unlist(JustSamples[[f]][names(JustSamples[[f]])== "file"]),nrow=length(unlist(JustSamples[[f]][names(JustSamples[[f]])== "file"]))/4,byrow=T)
                if(any(names(JustSamples[[f]]) == "file")){
                  if(any(grepl("Read 1 FASTQ",vecSample[,1]))){
                    vecSampleLocations[f] <- vecSample[grep("Read 1 FASTQ",vecSample[,1]),2]
                   if(any(grepl("Read 2 FASTQ",vecSample[,1]))){
                      vecSampleLocations[f] <- paste(vecSampleLocations[f],vecSample[grep("Read 2 FASTQ",vecSample[,1]),2],sep=";")
                    }
                    vecSampleNames[f] <- JustSamples[[f]]$name
                  }
                }


              }
              Temp2 <- cbind(matrix(rep(TempFull[j,],length(vecSampleNames)),ncol=ncol(TempFull),byrow=T),rep(Lanes[[l]]$lane,length(vecSampleNames)),vecSampleNames,vecSampleLocations)
                #Temp <- Temp[-1,]
              colnames(Temp2)  <- c("SLXID","Run","Type","Lane","SampleNames","Location")
              if(incount > 1){
                 Temp2Full <- rbind(Temp2Full,Temp2)
              }else{
                 Temp2Full <- Temp2
              }
           }
         }

       }

          if(j > 1){
            Temp3 <- rbind(Temp3,Temp2Full)
          }else{
            Temp3 <- Temp2Full
          }


    }
        if(i > 1){
          LimsLocations <- rbind(LimsLocations,Temp3)
        }else{
          LimsLocations <- Temp3
        }


  }
}

SLXIDLimsLocations <- LimsLocations

#SLXIDLimsLocations

write.table(SLXIDLimsLocations,LimsOutPut,sep="\t")
}

#findSLXIDS("SLX-7546","LimsInfo.txt")


