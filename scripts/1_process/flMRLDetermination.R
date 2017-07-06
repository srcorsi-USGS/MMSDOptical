# Script with workflow for determining minimum reporting levels (MRLs) for
# fluorescence for MMSD Phase III virus (Menomonee River)


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
Wavelength <- "exem" # column name defining excitation/emmission wavelengths

### Phase 4 data ###
# Load Phase 4 summary data, 3-D fluorescence, and absorbance data
load(file.path(raw.path,"PhaseIV","VirusPhaseIVData.Rdata"))
load(file.path(raw.path,"PhaseIV","MMSDOpticalData.RData"))
dfOptSum <- df
df <- read.csv(file.path(raw.path,"PhaseIV","MMSDOptSummary.csv"),skip=1,stringsAsFactors = FALSE)
names(dffl)[1] <- Wavelength

#truncate column names at 7 characters to match GR numbers.
names(dffl) <- substr(names(dffl),1,7)

#Remove rows with all NAs
naRows <- which(apply(dffl,1,FUN=function(x)mean(is.na(x)))>0.5)
dffl <- dffl[-naRows,]

#Define GRnumbers for environmental samples in filtered summary file
sampleGRnums <- dfOptSum[which(!is.na(dfOptSum[,GRnumber])),"GRnumber"]
sampleGRnums <- sampleGRnums[which(sampleGRnums %in% names(dffl))]

#Define GRnumbers for blank samples in original summary file from CA
blankRows <- grep("blank",df$QA, ignore.case = TRUE)
blankRows2 <- grep("blank",df$SampleNotes, ignore.case = TRUE)
sum(!(blankRows2 %in% blankRows)) #confirmed that all blanks are represented in blankRows
blankGRnums <- df[blankRows,GRnumber]
blankSummary <- df[blankRows,]
blankSummary <- subset(blankSummary,!is.na(A254))
blankGRnums <- blankGRnums[which(blankGRnums %in% names(dffl))]

dfblanks <- dffl[,c(Wavelength,blankGRnums)]
dfBlankLong <- gather(dfblanks,GRnumber,value,-exem)
blankOutliers <- unique(dfBlankLong[which(dfBlankLong$value > 0.04),"GRnumber"])

xnums <- as.numeric(factor(dfBlankLong$exem,levels=dffl$exem))

filenm <- "P4BlankPlots.pdf"
pdf(filenm,width = 20,height = 8)
for(exWave in unique(substr(dfblanks[,Wavelength],1,3))){

  #exWave <- 240
exemRange <- which(substr(dfBlankLong[,Wavelength],1,3)==exWave)
exemNames <- which(substr(dfblanks[,Wavelength],1,3)==exWave)
  subdf <- dfBlankLong[exemRange,]
  subdf <- subset(subdf,!is.na(value))
xnums <- as.numeric(factor(subdf$exem,levels=dfblanks[exemNames, "exem"]))
waveNames <- dfblanks[exemNames, "exem"]
plot(subdf$value~xnums,main=exWave,xaxt="n")
axis(side=1,at=1:length(waveNames),labels=waveNames,las=2,cex=0.5)

#subdf[which(subdf$value>0.08),Wavelength]
#subdf[which(subdf$value>0.12),"GRnumber"]
}
dev.off()
shell.exec(filenm)


which(dfblanks[which(dfblanks[,Wavelength]=="265/292"),]>0.2)

# Compute MRLs and adjust vectorized fl data
dfMRLs <- optMRL(dffl,Wavelength,blankGRnums)
flList <- optMRLAdjust(dffl,dfMRLs,Wavelength,sampleGRnums)

#DF of adjusted fl values
dffl2 <- flList[[1]]

saveRDS(dffl2,file=file.path(cached.path,cached.save,"flP4MRLAdjusted.rds"))
saveRDS(flList[[2]],file=file.path(cached.path,cached.save,"flP4withRemarks.rds"))
saveRDS(dfMRLs,file=file.path(cached.path,cached.save,"flMRLs.rds"))

# Explore results with a few graphics
plot(dfMRLs$Wavelength,dfMRLs$MRL)

MRLdiff <- dfMRLs[1:550,"MRL"] - dfMRLs[2:551,"MRL"]
plot(MRLdiff)
which(MRLdiff > 0.002)

#Examine the fraction of fluorescence values that are censored
#Most wavelengths that we are inerested in (file = ex_ems_means.csv in raw_data directory)
#are low or no percentage censored except peak B (275/304) and the wastewater range signal 
#which is a rectangle from ex = 255 to 290 to em = 302 to 350.

threshFraction <- 0.01
fractionCensored <- apply(flList[[2]][,-1],1,FUN=function(x) length(grep("<",x))/length(x))
names(fractionCensored) <- flList[[2]][,1]
wavesCensored <- fractionCensored[which(fractionCensored>threshFraction)]

filenm <- "fractionCensored.pdf"
pdf(filenm,width=15,height=8)
plot(wavesCensored,xaxt="n")
axis(side=1,at=1:length(wavesCensored),labels = names(wavesCensored),las=2,cex.axis=0.75)

dfwavesCensored <- separate(data=data.frame(exem=names(wavesCensored)),col = exem,into = c("ex","em"),sep = "/")
dfwavesCensored$censoredfraction <- wavesCensored
dfwavesCensored$ex <- as.numeric(dfwavesCensored$ex)
dfwavesCensored$em <- as.numeric(dfwavesCensored$em)

colfunc <- colorRampPalette(c("lightblue","darkblue"))
plotColors <- colfunc(10)
dfwavesCensored$plotColors <- plotColors[round(dfwavesCensored$censoredfraction*9)+1]

plot(dfwavesCensored$ex,dfwavesCensored$em,
     xlab="Excitation wavelength (nm)",ylab="Emmission wavelength(nm)",
     col=dfwavesCensored$plotColors,pch=20)
mtext(side=3,line=1,font=2,
  paste0("Excitation/Emmission Wavelengths with greater than a ",threshFraction," proportion censored values"))
legend("topleft",legend = paste("<",c(1:10)/10),col = plotColors,pch=20)

dev.off()
shell.exec(filenm)

# # Generate plots with original and adjusted values
# filenm <- "flPlotsP4.pdf"
# pdf(filenm)
# par(mfrow=c(2,1))
# for(flCol in sampleGRnums){
#   plot(dffl2[,Wavelength],dffl2[,flCol],
#        type="l",lty=1,col="blue",
#        xlab="Wavelength (nm)",
#        ylab="fluorescence",
#        main = flCol)
#   lines(dffl[,Wavelength],dffl[,flCol],
#         lty=3,col="orange")
# }
# dev.off()
# shell.exec(filenm)
# 
# 
# #Generate plots of the difference between adjusted and original values
# filenm <- "flDiffPlotsP4.pdf"
# pdf(filenm)
# par(mfrow=c(2,1))
# for(flCol in sampleGRnums){
#   flDiff <- dffl2[,flCol] - dffl[,flCol]
#   plot(dffl2[,Wavelength],flDiff,
#        type="l",lty=1,col="blue",
#        xlab="Wavelength (nm)",
#        ylab="Fluorescence difference",
#        main = flCol)
# }
# dev.off()
# shell.exec(filenm)


### Phase 3 data ###
# Load Phase 3 summary data, 3-D fluorescence, and absorbance data

load(file.path(raw.path,"PhaseIII","MMSDabsEEMs.RData"))
load(file.path(raw.path,"PhaseIII","dfOptAnalysisDataMMSDJan2015.RData"))
dfMRLs <- readRDS(file.path(cached.path,cached.save,"flMRLs.rds"))
dffl <- dfFluor
Wavelength <- "exem"
names(dffl)[1] <- Wavelength

#Remove rows with all NAs
naRows <- which(apply(dffl,1,FUN=function(x)mean(is.na(x)))>0.5)
dffl <- dffl[-naRows,]
dfOptSum <- dfOptSumAll

#Define GRnumbers for environmental samples in filtered summary file
sampleGRnums <- dfOptSum[which(!is.na(dfOptSum[,GRnumber])),"GRnumber"]
sampleGRnums %in% names(dffl)

# Adjust vectorized fl data
flList <- optMRLAdjust(dffl,dfMRLs,Wavelength,sampleGRnums)

#DF of adjusted fl values
dffl2 <- flList[[1]]

saveRDS(dffl2,file=file.path(cached.path,cached.save,"flP3MRLAdjusted.rds"))
saveRDS(flList[[2]],file=file.path(cached.path,cached.save,"flP3withRemarks.rds"))

#Examine the fraction of fluorescence values that are censored
#Most wavelengths that we are inerested in (file = ex_ems_means.csv in raw_data directory)
#are low or no percentage censored except peak B (275/304) and the wastewater range signal 
#which is a rectangle from ex = 255 to 290 to em = 302 to 350.

threshFraction <- 0.01
fractionCensored <- apply(flList[[2]][,-1],1,FUN=function(x) length(grep("<",x))/length(x))
names(fractionCensored) <- flList[[2]][,1]
wavesCensored <- fractionCensored[which(fractionCensored>threshFraction)]

filenm <- "fractionCensoredP3.pdf"
pdf(filenm,width=15,height=8)
plot(wavesCensored,xaxt="n")
axis(side=1,at=1:length(wavesCensored),labels = names(wavesCensored),las=2,cex.axis=0.75)

dfwavesCensored <- separate(data=data.frame(exem=names(wavesCensored)),col = exem,into = c("ex","em"),sep = "/")
dfwavesCensored$censoredfraction <- wavesCensored
dfwavesCensored$ex <- as.numeric(dfwavesCensored$ex)
dfwavesCensored$em <- as.numeric(dfwavesCensored$em)

colfunc <- colorRampPalette(c("lightblue","darkblue"))
plotColors <- colfunc(10)
dfwavesCensored$plotColors <- plotColors[round(dfwavesCensored$censoredfraction*9)+1]

plot(dfwavesCensored$ex,dfwavesCensored$em,
     xlab="Excitation wavelength (nm)",ylab="Emmission wavelength(nm)",
     col=dfwavesCensored$plotColors,pch=20)
mtext(side=3,line=1,font=2,
      paste0("Excitation/Emmission Wavelengths with greater than a ",threshFraction," proportion censored values"))
legend("topleft",legend = paste("<",c(1:10)/10),col = plotColors,pch=20)

dev.off()
shell.exec(filenm)


# 
# # Generate plots with original and adjusted values
# filenm <- "flPlotsP3.pdf"
# pdf(filenm)
# par(mfrow=c(2,1))
# for(flCol in sampleGRnums){
#   plot(dffl2[,Wavelength],dffl2[,flCol],
#        type="l",lty=1,col="blue",
#        xlab="Wavelength (nm)",
#        ylab="Fluorescence",
#        main = flCol)
#   lines(dffl[,Wavelength],dffl[,flCol],
#         lty=3,col="orange")
# }
# dev.off()
# shell.exec(filenm)
# 
# 
# #Generate plots of the difference between adjusted and original values
# filenm <- "flDiffPlotsP3.pdf"
# pdf(filenm)
# par(mfrow=c(2,1))
# for(flCol in sampleGRnums){
#   flDiff <- dffl2[,flCol] - dffl[,flCol]
#   plot(dffl2[,Wavelength],flDiff,
#        type="l",lty=1,col="blue",
#        xlab="Wavelength (nm)",
#        ylab="Fluorescence difference",
#        main = flCol)
# }
# dev.off()
# shell.exec(filenm)
