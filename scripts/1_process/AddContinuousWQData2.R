# Gather information from continuous water quality sensors for matching with
# optical and bacteria data

library(dplyr)
library(USGSHydroTools)

#set data directories
munge.path <- "0_munge"
WQ.path <- "ContinuousWQ"
cached.path <- "cached_data"
summary.path <- "SummaryVariables"
summary.save <- "1_SummaryVariables"
cached.save <- "0_munge"

#Read dataframe with bacteria and optical data
df <- readRDS(file.path(cached.path,summary.save,"dfOptP3P4Combined.rds"))
df.orig <- df

sites <- c("04087030","04087050","04087088","04087119","04087120","04087142","04086600","05426060")
names(sites) <- c("MF","LD","UW","HW","MW","MC","CG","BK")
# parms <- c("00060","00095","00300","00010","63680")
# names(parms) <- c("Q","SC","DO","WT","Turbidity")
parmNames <- c("Discharge","Dissolved_oxygen","Specific_cond","Turbidity","Water_Temperature")

dfSum <- df

#Add summary stats for each continuous water quality varible
for(j in c(1,2,3,4,5)){
  parm <- parmNames[j]
  
  #Read UV data file
  filenm <- paste0(parm,"AllSites.rds")
  dfWQ <- readRDS(file.path(cached.path,munge.path, WQ.path,filenm))
  dfList <- list()
  for(i in 1:6){
    site <- names(sites[i])
    subdf <- subset(dfWQ,STAID == sites[i])
    subdfSum <- subset(dfSum,abbrev==names(sites[i]))

    #Add columns that represent the stats for each individual WQ parm
    dfList[[i]] <- TSstormstats(df = dfWQ,date="pdate",varname="Value",dates=subdfSum,
                                starttime="psdate",endtime="pedate",
                                stats.return = c("mean","median","max","min"),
                                out.varname = parm)
  }
  
  #Combine individual sites back into one file
  dfSum <- dfList[[1]]
  for(k in 2:length(dfList)){
    dfSum <- rbind(dfSum,dfList[[k]])
  }
}

saveRDS(dfSum,file.path(cached.path,summary.save,"dfOptP3P4CombinedContWQ.rds"))
