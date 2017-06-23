# Script with workflow for adding summary optical data from vectorized
# fluorescence and full absorbance scans

###################
#Set working directory
SNsys <- system("wmic bios get serialnumber", intern = TRUE, show.output.on.console = FALSE)
SN <- gsub("\\s","",SNsys)[2]
baseDir <- "MMSD/Phase IV/virus/Data/Optical"
if(SN == "PB519R1") {Project <- paste("D:/srcldata/",baseDir,sep="")
}else {Project <- paste("M:/QW Monitoring Team/",baseDir,sep="")
}

setwd(Project)
###################

library(USGSHydroOpt)
library(data.table)

checkDups <- function(df,parm){
  df[duplicated(df[,parm]),parm]
}

# Load fluorescence and absorbnce data

##########################################################################
# dffl <- fread('MMSDFlmx4.csv',sep = ',',stringsAsFactors = FALSE)
# dfabs <- fread('MMSDAbs.csv')
# ## Generate 3d array of EEMs data ##
# MMSDP43DEEMs <-VectorizedTo3DArray(df = dffl,ExEm="ex/em",grnum='GRnumber')
# save(dffl,dfabs,MMSDP43DEEMs,file='MMSDOpticalData.RData')
# 
load('MMSDOpticalData.RData')
# dfOpt <- read.csv('MMSDOptSummary.csv',stringsAsFactors = FALSE,skip=1)
# 
# dfOpt1 <- dfOpt[,-(which(names(dfOpt)=="A254"):dim(dfOpt)[2])]
load(file='../merged/VirusPhaseIVData.Rdata')
dfOpt1 <- df

dfOpt1$GRnumber %in% substr(names(dffl),1,7)

dfCheck <- dfOpt1[dfOpt1$GRnumber %in% names(dffl),]

#dfOpt1$GRnumber %in% names(dfabs)
load("M:/QW Monitoring Team/GLRI toxics/Shared_optical_data/GLRI/Final data/ratioOrder2016-03-03.Rdata") # Ratio orders determined by GLRI data set

# Read summary signals to extract from Fl and abs info
#SummaryDir <- "D:/SRCLData/AqualogProcessing/SummaryVariables/MMSD/"
SummaryDir <- "M:/QW Monitoring Team/Optical sensors/AqualogProcessing/SummaryVariables/MMSD/"

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


# Add summary variables to summary data frame

##############################################################################################
######### Add summary variables to summary data frame #######################################

#Fluorescence pairs and means
dfOpt2 <- getMeanFl(a=MMSDP43DEEMs,signals=dfFlSignals,Peak="Peak",Ex1="Ex1",Ex2="Ex2",Em1="Em1",Em2="Em2",dataSummary=dfOpt1,grnum="GRnumber")

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


fDate <- as.character(Sys.Date())
save(dfOptSumAll,file=paste("dfOptSummaryMMSDP3",fDate,".RData",sep=""))
write.csv(dfOptSumAll,paste("dfOptSummaryMMSDP3",fDate,".csv",sep=""),row.names=FALSE)

varNames <- names(dfOptSumAll)
varNames
write.csv(names(dfOptSumAll),file=paste("namesAllMMSDP3",fDate,".csv",sep=""),row.names=FALSE)




