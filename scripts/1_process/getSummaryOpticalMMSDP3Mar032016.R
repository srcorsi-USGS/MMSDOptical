# Script with workflow for adding summary optical data from vectorized
# fluorescence and full absorbance scans for MMSD Phase III virus (Menomonee River)

library(USGSHydroOpt)

#set data directories
raw.path <- "raw_data"
cached.path <- "cached_data"
summary.path <- "SummaryVariables"
summary.save <- "1_SummaryVariables"
cached.save <- "0_munge"

checkDups <- function(df,parm){
  df[duplicated(df[,parm]),parm]
}

# Load summary data, 3-D fluorescence, and absorbance data

load(file.path(raw.path,"PhaseIII","MMSDabsEEMs.RData"))
load(file.path(raw.path,"PhaseIII","MMSD3DEEMs.RData"))
load(file.path(raw.path,"PhaseIII","dfOptAnalysisDataMMSDJan2015.RData"))
dfOpt1 <- dfOptSumAll[,-(which(names(dfOptSumAll)=="OB1"):dim(dfOptSumAll)[2])]
#dfOpt1$GRnumber %in% names(dfabs)

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
dfOpt2 <- getMeanFl(a=MMSD3DEEMs,signals=dfFlSignals,Peak="Peak",Ex1="Ex1",Ex2="Ex2",Em1="Em1",Em2="Em2",dataSummary=dfOpt1,grnum="GRnumber")

#HIX, FI, Freshness 
dfOpt2 <- getIndexes(a=MMSD3DEEMs,dataSummary=dfOpt2,grnum="GRnumber")

#Single absorbance signals
dfOpt2 <- getAbs(dataAbs=dfabs,waveCol="Wavelength",wavs=dfAbsSignals[,1],colSubsetString="gr",dataSummary=dfOpt2,grnum="GRnumber")

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

saveRDS(dfOptSumAll,file=file.path(cached.path,summary.save,"dfOptSummaryMMSDP3.rds"))





