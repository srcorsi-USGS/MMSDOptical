# Script with workflow for adding summary optical data from vectorized
# fluorescence and full absorbance scans for MMSD Phase IV virus (Menomonee, Milwaukee, Bark Rivers)

library(USGSHydroOpt)
library(data.table)
library(tidyr)

#set data directories
raw.path <- "raw_data"
cached.path <- "cached_data"
munge.path <- "0_munge"
summary.path <- "SummaryVariables"
summary.save <- "1_SummaryVariables"
cached.save <- "0_munge"
processed.path <- "0_ProcessMDLs"

checkDups <- function(df,parm){
  df[duplicated(df[,parm]),parm]
}

# Load fluorescence and absorbnce data

##########################################################################
# Load summary data, 3-D fluorescence, and absorbance data
# Load summary data, vectorized fluorescence, and absorbance data
df <- readRDS(file.path(cached.path,munge.path,"VirusPhaseIVData.rds"))
dfabs <- readRDS(file.path(cached.path,processed.path, "absP4MRLAdjusted.rds"))
dffl <- readRDS(file.path(cached.path,processed.path, "flP4MRLAdjusted.rds"))

#Generate 3-D EEMS array
MMSDP4DEEMs <- VectorizedTo3DArray(dffl,"exem", "GRnumber")

#Name the wavelength column in the abs file
names(dfabs)[1] <- "Wavelength"

dfOpt1 <- df[,-(which(names(df)=="A254"):dim(df)[2])]
#dfOpt1$GRnumber %in% names(dfabs)

which(!(dfOpt1$GRnumber %in% names(dfabs)))
dfOpt1[which(!(dfOpt1$GRnumber %in% substr(names(dffl),1,7))),"GRnumber"]

#Missing samples in here are included in the Aqualog analysis but not the fluoromax
dfMissingFl <- dfOpt1[which(!(dfOpt1$GRnumber %in% substr(names(dffl),1,7))),]

#Subset the summary file to include only those that are present in the fluorescence file
dfOpt1 <- dfOpt1[which((dfOpt1$GRnumber %in% substr(names(dffl),1,7))),]


names(dffl) <- substr(names(dffl),1,7)

# Read summary signals to extract from Fl and abs info
SummaryDir <- paste0("./",raw.path,"/",summary.path,"/")

dfFlSignals <- read.csv(paste(SummaryDir,"ex_ems_means.csv",sep=""),as.is=TRUE)
dfAbsSignals <- read.csv(paste(SummaryDir,"abs_wavs.csv",sep=""),as.is=TRUE)
AbsSignals <- as.numeric(dfAbsSignals[,1])
dfSagSignals <- read.csv(paste(SummaryDir,"SagWaves.csv",sep=""),as.is=TRUE)
ratioSignalsAbs <- read.csv(paste(SummaryDir,"ratioSignalsAbs.csv",sep=""),as.is=TRUE)
ratioSignalsAbs <- ratioSignalsAbs[which(ratioSignalsAbs[2]>0),1]
ratioSignalsSr <- read.csv(paste(SummaryDir,"ratioSignalsSrCA.csv",sep=""),as.is=TRUE)
ratioSignalsSr <- ratioSignalsSr[which(ratioSignalsSr[2]>0),1]
ratioSignalsSniff <- read.csv(paste(SummaryDir,"ratioSignalsSniff.csv",sep=""),as.is=TRUE)
ratioSignalsSniff <- ratioSignalsSniff[which(ratioSignalsSniff[2]>0),1]
logSignals <- read.csv(paste(SummaryDir,"logSignals.csv",sep=""),as.is=TRUE)[,1]
ratioOrder <- readRDS(file.path(cached.path,summary.save,"ratioOrder.rds"))


# Add summary variables to summary data frame

##############################################################################################
######### Add summary variables to summary data frame #######################################

#Fluorescence pairs and means
dfOpt2 <- getMeanFl(a=MMSDP4DEEMs,signals=dfFlSignals,Peak="Peak",Ex1="Ex1",Ex2="Ex2",Em1="Em1",Em2="Em2",dataSummary=dfOpt1,grnum="GRnumber")

#HIX, FI, Freshness 
dfOpt2 <- getIndexes(a=MMSDP4DEEMs,dataSummary=dfOpt2,grnum="GRnumber")

#Single absorbance signals
dfOpt2 <- getAbs(dataAbs=dfabs,waveCol="Wavelength",wavs=dfAbsSignals[,1],colSubsetString="gr",dataSummary=dfOpt2,grnum="GRnumber")

#START HERE....PROBLEM WITH MISSING VALUE ON 
#Spectral slopes
dfOpt2 <- getSag(dataAbs=dfabs,waveCol="Wavelength",sag=dfSagSignals,colSubsetString="gr",dataSummary=dfOpt2,grnum="GRnumber")

#deviance of abs from exponential regression in the wastewater area
dfOpt2 <- getExpResid(wavelength=267,rangeReg=c(240,340),rangeGap=c(255,300),dataAbs=dfabs,waveCol="Wavelength",colSubsetString="gr",dataSummary=dfOpt2,grnum="GRnumber")

#write.csv(names(dfOpt2),file="dfOptNames.csv")

#Ratios of a few things
dfOpt2 <- getRatios(dataSummary=dfOpt2,grnum="GRnumber",specifyOrder = TRUE,ratioVars = ratioOrder)

#log transform where it makes sense
dfOpt2 <- getLog10(dataSummary=dfOpt2,signals=logSignals,grnum="GRnumber")

##############################################################################################
##############################################################################################

dfOptSumAll <- dfOpt2 #merge(dfAll,dfOpt2,by=c("GRnumber","Site"),)
checkDups(dfOptSumAll,"GRnumber")
#checkDups(dfAll,"GRnumber")

#Check for potential issues
nonFinite <- function(x) {sum(!is.finite(x))}

bad <- numeric()
for(var in names(dfOpt2)){
  bad <- c(bad,nonFinite(dfOpt2[,var]))
}
names(bad) <- names(dfOpt2)
bad[which(bad>0)]

# Determine ratio order and save it for use with Phase III summary optical determination
names(dfOptSumAll)
ratioRows <- which(substr(names(dfOptSumAll),start = 1,stop = 1)=="r")
ratioVars <- substr(names(dfOptSumAll)[ratioRows],2,50)

r1 <- character()
r2 <- character()
for(i in 1:length(ratioVars)){
  if(substr(ratioVars[i],1,3) == "Sag"){
    r1 <- c(r1,substr(ratioVars[i],1,10))
    r2 <- c(r2,substr(ratioVars[i],12,50))
  }else{
    r1 <- c(r1,strsplit(ratioVars[i],"_")[[1]][1])
    r2 <- c(r2,strsplit(ratioVars[i],"_")[[1]][2])
  }
}

dfRatioOrder <- data.frame(var1 = r1,var2 = r2,stringsAsFactors = FALSE)

saveRDS(dfOptSumAll,file=file.path(cached.path,summary.save,"dfOptSummaryMMSDP4.rds"))
        
