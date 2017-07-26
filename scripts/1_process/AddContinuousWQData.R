# Gather information from continuous water quality sensors for matching with
# optical and bacteria data

library(dataRetrieval)
library(dplyr)

#set data directories
raw.path <- "raw_data"
WQ.path <- "ContinuousWQ"
cached.path <- "cached_data"
summary.path <- "SummaryVariables"
summary.save <- "1_SummaryVariables"
cached.save <- "0_munge"


df <- readRDS(file.path(cached.path,summary.save,"dfOptP3P4Combined.rds"))
df.orig <- df

#Define sites and parameters
sites <- c("04087030","04087050","04087088","04087119","04087120","04087142","04086600","05426060")
names(sites) <- c("MF","LD","UW","HW","MW","MC","CG","BK")
parms <- c("00060","00095","00300","00010","63680")
names(parms) <- c("Q","SC","DO","WT","Turbidity")

parmNames <- c("Discharge","Dissolved_oxygen","Specific_cond","Turbidity","Water_Temperature")

for(j in 1: length(parmNames)){
  filenm <- grep(parmNames[j],list.files(path = file.path(raw.path,WQ.path,parmNames[j])),value=TRUE)
  for(i in 1:6){
    fileNum <- grep(sites[i],filenm)
    dfWQsub <- read.csv(
      file.path(raw.path,WQ.path,parmNames[j],filenm[fileNum]),
      skip=14,stringsAsFactors = FALSE)
    dfWQsub$STAID <- sites[i]
    
    if(i == 1){dfWQ <- dfWQsub
    }else dfWQ <- rbind(dfWQ,dfWQsub)
  }
  
  dfWQ$pdate <- as.POSIXct(dfWQ$Timestamp..UTC.06.00.,format="%Y-%m-%d %H:%M:%S", tz="CST6CDT")
  dfWQ$pdate <- as.POSIXct(format(as.POSIXct(dfWQ$pdate),tz="Etc/GMT+6",usetz=TRUE),tz="Etc/GMT+6")
  saveRDS(dfWQ,file=file.path(cached.path,cached.save,WQ.path,paste0(parmNames[j],"AllSites.rds")))
}



