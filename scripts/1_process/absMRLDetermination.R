# Script with workflow for adding summary optical data from vectorized
# fluorescence and full absorbance scans for MMSD Phase III virus (Menomonee River)

library(USGSHydroOpt)

#set data directories
raw.path <- "raw_data"
cached.path <- "cached_data"
summary.path <- "SummaryVariables"
summary.save <- "1_SummaryVariables"
cached.save <- "0_munge"

GRnumber <- "GRnumber"
# Load summary data, 3-D fluorescence, and absorbance data

# load(file.path(raw.path,"PhaseIII","MMSDabsEEMs.RData"))
# load(file.path(raw.path,"PhaseIII","MMSD3DEEMs.RData"))
# load(file.path(raw.path,"PhaseIII","dfOptAnalysisDataMMSDJan2015.RData"))

load(file.path(raw.path,"PhaseIV","VirusPhaseIVData.Rdata"))
load(file.path(raw.path,"PhaseIV","MMSDOpticalData.RData"))
dfOptSum <- df
df <- read.csv(file.path(raw.path,"PhaseIV","MMSDOptSummary.csv"),skip=1,stringsAsFactors = FALSE)
names(dfabs)[1] <- "Wavelength"

#Define GRnumbers for environmental samples in filtered summary file
sampleGRnums <- dfOptSum[which(!is.na(dfOptSum$GRnumber)),"GRnumber"]
sampleCols <- which(names(dfabs) %in% sampleGRnums)

#Define blank samples in original summary file from CA
blankRows <- grep("blank",df$QA, ignore.case = TRUE)
blankRows2 <- grep("blank",df$SampleNotes, ignore.case = TRUE)
sum(!(blankRows2 %in% blankRows)) #confirmed that all blanks are represented in blankRows
blankGRnums <- df[blankRows,GRnumber]
blankCols <- which(names(dfabs) %in% blankGRnums)

# Generate data frame with information on the blank samples by wavelength and
# compute the minimum reporting level based on mean + 3 * SD for the blank samples
# If the mean is less than zero, set the MRL to 3* SD
dfBlankSummary <- data.frame(Wavelength = dfabs$Wavelength)
dfBlankSummary$mean <- apply(dfabs[,blankGRnums],MARGIN = 1, mean, na.rm=TRUE)
dfBlankSummary$min <- apply(dfabs[,blankGRnums],MARGIN = 1, min, na.rm=TRUE)
dfBlankSummary$sd <- apply(dfabs[,blankGRnums],MARGIN = 1, sd, na.rm=TRUE)
dfBlankSummary$MRL <- dfBlankSummary$mean + 3 * dfBlankSummary$sd
dfBlankSummary$MRL <- ifelse(dfBlankSummary$mean < 0, 3 * dfBlankSummary$sd, dfBlankSummary$MRL)

#Generate data frame with adjusted values based on the MRL. Generate a second 
#dataframe with remark columns indicating values that are less than the MRL
dfabs2 <- data.frame(Wavelength = dfabs[,"Wavelength"])
dfabsRemarks <- data.frame(Wavelength = dfabs[,"Wavelength"])
for(i in sampleCols){
  dfabs2 <- cbind(dfabs2,ifelse(dfabs[,i] < dfBlankSummary[,"MRL"],dfBlankSummary[,"MRL"],dfabs[,i]))
  dfabsRemarks <- cbind(dfabsRemarks,ifelse(dfabs[,i] < dfBlankSummary[,"MRL"],paste("<",dfBlankSummary[,"MRL"]),dfabs[,i]))
}

names(dfabs2)[-1] <- names(dfabs)[sampleCols]
names(dfabsRemarks)[-1] <- names(dfabs)[sampleCols]

filenm <- "absPlots.pdf"
pdf(filenm)
par(mfrow=c(2,1))
for(absCol in sampleGRnums){
plot(dfabs2[,"Wavelength"],dfabs2[,absCol],
     type="l",lty=1,col="blue",
     xlab="Wavelength (nm)",
     ylab="Absorbance coefficient",
     main = absCol)
  lines(dfabs[,"Wavelength"],dfabs[,absCol],
       lty=3,col="orange")
}
dev.off()
shell.exec(filenm)

filenm <- "absDiffPlots.pdf"
pdf(filenm)
par(mfrow=c(2,1))
for(absCol in sampleGRnums){
  absDiff <- dfabs2[,absCol] - dfabs[,absCol]
  plot(dfabs2[,"Wavelength"],absDiff,
       type="l",lty=1,col="blue",
       xlab="Wavelength (nm)",
       ylab="Absorbance difference",
       main = absCol)
}
dev.off()
shell.exec(filenm)

