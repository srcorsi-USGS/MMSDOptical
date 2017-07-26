# Gather information from continuous water quality sensors for matching with
# optical and bacteria data

library(dataRetrieval)

#set data directories
raw.path <- "raw_data"
WQ.path <- "ContinuousWQ"
cached.path <- "cached_data"
summary.path <- "SummaryVariables"
summary.save <- "1_SummaryVariables"
cached.save <- "0_munge"

#Define function to query multiple sites for UV data. Return a list of dataframes,
# one per site.
getUVList <- function(sites,pcodes,startDate,endDate,tz="UTC"){
  UVlist <- list()
  for(i in 1:length(sites)){
    setAccess("internal")
    site <- sites[i]
    UVlist[[i]] <- readNWISuv(siteNumbers = site,parameterCd = pcode, startDate = startDate,endDate = endDate,tz=tz)
  }
  return(UVlist)
}


df <- readRDS(file.path(cached.path,summary.save,"dfOptP3P4Combined.rds"))
df.orig <- df



dfTurb <- read.csv(
  file.path(raw.path,WQ.path,"Turbidity,_Form_Neph._FNU.YSI.Work@04087120.20101001.csv"),
  skip=14)
dfTurb$pdate <- as.POSIXct(dfTurb$Timestamp..UTC.06.00.,format="%Y-%m-%d %H:%M:%S", tz="CST6CDT")
dfTurb$pdate <- as.POSIXct(dfTurb$pdate, tz="Etc/GMT-6")
"Discharge.ft^3_s@04087120.20101001.csv"
sites <- c("04087030","04087050","04087088","04087119","04087120","04087142","04086600","05426060")
names(sites) <- c("MF","LD","UW","HW","MW","MC","CG","BK")
parms <- c("00060","00095","00300","00010","63680")
names(parms) <- c("Q","SC","DO","WT","Turbidity")

startDate <-"2011-01-01" 
endDate <- "2014-07-30"
tz <- "UTC"

# #Retrieve Q for all sites
# i <- 1
# Qlist <- getUVList(sites = sites,pcode = parms[i],startDate = startDate,endDate = endDate,tz = tz)
# names(Qlist) <- names(sites)
# saveRDS(Qlist,file.path(raw.path,"Qlist_2010_2015.rds"))
        
#Retrieve SC for all sites
i <- 2
SClist <- getUVList(sites,parms[i],startDate,endDate,tz)
Qlistnames(SClist) <- names(sites)
saveRDS(SClist,file.path(raw.path,"Qlist_2010_2015.rds"))

range(df$pedate)

