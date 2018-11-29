# Script with workflow for determining minimum reporting levels (MRLs) for
# absorbance for MMSD Phase III virus (Menomonee River)


# Workflow for determing MRLs and adjusting data for Phase III and IV data

library(USGSHydroOpt)
library(tidyr)

#set data directories
raw.path <- "raw_data"
cached.path <- "cached_data"
summary.path <- "SummaryVariables"
summary.save <- "1_SummaryVariables"
cached.save <- "0_ProcessMDLs"
script.path <- "scripts"
process.path <- "1_process"

#Load functions for determining MRLs and adjusting data based on MRLs
source(file.path(script.path,process.path,"optMRL.R"))
source(file.path(script.path,process.path,"optMRLAdjust.R"))

GRnumber <- "GRnumber" #Lab ID number
Wavelength <- "Wavelength" # column name defining absorbance wavelengths

### Phase 4 data ###
# Load Phase 4 summary data, 3-D fluorescence, and absorbance data
load(file.path(raw.path,"PhaseIV","VirusPhaseIVData.Rdata"))
load(file.path(raw.path,"PhaseIV","MMSDOpticalData.RData"))
dfOptSum <- df
df <- read.csv(file.path(raw.path,"PhaseIV","MMSDOptSummary.csv"),skip=1,stringsAsFactors = FALSE)
names(dfabs)[1] <- Wavelength
dfabs <- as.data.frame(dfabs)

#Define GRnumbers for environmental samples in filtered summary file
sampleGRnums <- dfOptSum[which(!is.na(dfOptSum[,GRnumber])),"GRnumber"]

#Define GRnumbers for blank samples in original summary file from CA
blankRows <- grep("blank",df$QA, ignore.case = TRUE)
blankRows2 <- grep("blank",df$SampleNotes, ignore.case = TRUE)
sum(!(blankRows2 %in% blankRows)) #confirmed that all blanks are represented in blankRows
blankGRnums <- df[blankRows,GRnumber]

dfblanks <- dfabs[,c(Wavelength,blankGRnums)]
dfBlankLong <- gather(dfblanks,GRnumber,value,-Wavelength)
blankOutliers <- unique(dfBlankLong[which(dfBlankLong$value > 0.04),"GRnumber"])
plot(dfBlankLong$value~as.numeric(dfBlankLong$Wavelength),xlim=c(250,300),ylim=c(-0.01,0.01))


# Compute MRLs and adjust vectorized abs data
dfMRLs <- optMRL(dfabs,Wavelength,blankGRnums)
absList <- optMRLAdjust(dfabs,dfMRLs,Wavelength,sampleGRnums)

#DF of adjusted abs values
dfabs2 <- absList[[1]]

dir.create(file.path(cached.path,cached.save), showWarnings = FALSE)

saveRDS(dfabs2,file=file.path(cached.path,cached.save,"absP4MRLAdjusted.rds"))
saveRDS(absList[[2]],file=file.path(cached.path,cached.save,"absP4withRemarks.rds"))
saveRDS(dfMRLs,file=file.path(cached.path,cached.save,"absMRLs.rds"))

# Explore results with a few graphics
plot(dfMRLs$Wavelength,dfMRLs$MRL)

MRLdiff <- dfMRLs[1:550,"MRL"] - dfMRLs[2:551,"MRL"]
plot(MRLdiff)
which(MRLdiff > 0.002)

# Generate plots with original and adjusted values
filenm <- "absPlotsP4.pdf"
pdf(filenm)
par(mfrow=c(2,1))
for(absCol in sampleGRnums){
  plot(dfabs2[,Wavelength],dfabs2[,absCol],
       type="l",lty=1,col="blue",
       xlab="Wavelength (nm)",
       ylab="Absorbance coefficient",
       main = absCol)
  lines(dfabs[,Wavelength],dfabs[,absCol],
        lty=3,col="orange")
}
dev.off()
#shell.exec(filenm)


#Generate plots of the difference between adjusted and original values
filenm <- "absDiffPlotsP4.pdf"
pdf(filenm)
par(mfrow=c(2,1))
for(absCol in sampleGRnums){
  absDiff <- dfabs2[,absCol] - dfabs[,absCol]
  plot(dfabs2[,Wavelength],absDiff,
       type="l",lty=1,col="blue",
       xlab="Wavelength (nm)",
       ylab="Absorbance difference",
       main = absCol)
}
dev.off()
#shell.exec(filenm)


### Phase 3 data ###
# Load Phase 3 summary data, 3-D fluorescence, and absorbance data

load(file.path(raw.path,"PhaseIII","MMSDabsEEMs.RData"))
load(file.path(raw.path,"PhaseIII","dfOptAnalysisDataMMSDJan2015.RData"))
dfMRLs <- readRDS(file.path(cached.path,cached.save,"absMRLs.rds"))

dfOptSum <- dfOptSumAll

#Define GRnumbers for environmental samples in filtered summary file
sampleGRnums <- dfOptSum[which(!is.na(dfOptSum[,GRnumber])),"GRnumber"]

# Adjust vectorized abs data
absList <- optMRLAdjust(dfabs,dfMRLs,Wavelength,sampleGRnums)

#DF of adjusted abs values
dfabs2 <- absList[[1]]

saveRDS(dfabs2,file=file.path(cached.path,cached.save,"absP3MRLAdjusted.rds"))
saveRDS(absList[[2]],file=file.path(cached.path,cached.save,"absP3withRemarks.rds"))


# Generate plots with original and adjusted values
filenm <- "absPlotsP3.pdf"
pdf(filenm)
par(mfrow=c(2,1))
for(absCol in sampleGRnums){
  plot(dfabs2[,Wavelength],dfabs2[,absCol],
       type="l",lty=1,col="blue",
       xlab="Wavelength (nm)",
       ylab="Absorbance coefficient",
       main = absCol)
  lines(dfabs[,Wavelength],dfabs[,absCol],
        lty=3,col="orange")
}
dev.off()
#shell.exec(filenm)


#Generate plots of the difference between adjusted and original values
filenm <- "absDiffPlotsP3.pdf"
pdf(filenm)
par(mfrow=c(2,1))
for(absCol in sampleGRnums){
  absDiff <- dfabs2[,absCol] - dfabs[,absCol]
  plot(dfabs2[,Wavelength],absDiff,
       type="l",lty=1,col="blue",
       xlab="Wavelength (nm)",
       ylab="Absorbance difference",
       main = absCol)
}
dev.off()
#shell.exec(filenm)