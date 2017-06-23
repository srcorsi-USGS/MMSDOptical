# MMSD Phase IV virus sampling
# Merge bacteria and virus with sample tracking information
# Compute pathogen category summations
# Process bacteria BLD samples

#setwd ("M:/QW Monitoring Team/MMSD/Phase IV/virus/R")
library(car)
library(ggplot2)
source('fxn_multiGrep.R')
sumNAs <- function(x) sum(is.na(x))

#set data directories
raw.path <- "raw_data"
cached.path <- "cached_data"
cached.save <- "1_munge"

########################
#UWM Bacteria data
dfBact <- read.csv(file.path(raw.path,"PhaseIV","USGS_MMSD results-7-16-15.csv"),stringsAsFactors = FALSE)
LODs <- c(225,225,225,225,225)
bactAbbrev <- c('lachno','bacHum','ent','eColi','fc')
names(LODs) <- bactAbbrev
bacteriaNames <- c("Lachno_2.CN.100.ml","Bac.human.CN.100.ml",
                   "Enterococcus.CN.100.ml","E.coli.CN.100.ml" ,
                   "Fecal.coliforms.CFU.100.ml")
for(i in 1:length(bacteriaNames)){
  dfBact[,bactAbbrev[i]] <- ifelse(dfBact[,bacteriaNames[i]]=="0",LODs[i],dfBact[,bacteriaNames[i]])
  dfBact[,bactAbbrev[i]] <- ifelse(dfBact[,bactAbbrev[i]]=="BLD",LODs[i],dfBact[,bactAbbrev[i]])
  dfBact[,bactAbbrev[i]] <- as.numeric(ifelse(dfBact[,bactAbbrev[i]]=="No DNA",LODs[i],dfBact[,bactAbbrev[i]]))
  
}

apply(dfBact[,bacteriaNames],2,sumNAs)


#form <- formula(paste('~',paste(bactAbbrev,collapse="+")))
#scatterplotMatrix(form, data=dfBact,main="Bacteria comparisons",log='xy')

########################
# USDA Virus data
dfVirus <- read.csv(file.path(raw.path,"PhaseIV","2012-2014 MMSD dataSRCRevised.csv"),stringsAsFactors = FALSE)
humanViruses <- c("Adenovirus.C.D.F","Adenovirus.A","Adenovirus.B","Enterovirus",
                  "G1.Norovirus","G2.Norovirus","Polyomavirus")
bovineViruses.orig <- c("Rotavirus.A","Coronavirus","Enterovirus.1","Adenovirus","Polyomavirus.1")
bovineViruses <- c("BRotavirus.A","BCoronavirus","BEnterovirus","BAdenovirus","BPolyomavirus")
Bacteria <- c("EHEC.eae.gene","EHEC.stx1.gene","EHEC.stx2.gene",
              "Salmonella.invA.gene","Salmonella.ttr.gene","Campylobacter.jejuni")
allPathogens <- c(humanViruses,bovineViruses,Bacteria)

names(dfVirus)[match(c(humanViruses,bovineViruses.orig,Bacteria),names(dfVirus))] <- 
  c(humanViruses,bovineViruses,Bacteria)
dfVirus$humanVirus <- apply(dfVirus[,humanViruses],1,sum)
dfVirus$bovineVirus <- apply(dfVirus[,bovineViruses],1,sum)
dfVirus$bacteria <- apply(dfVirus[,Bacteria],1,sum)
dfVirus$human <- apply(dfVirus[,c(humanViruses,Bacteria)],1,sum)
dfVirus$pathogens <- apply(dfVirus[,c(humanViruses,bovineViruses,Bacteria)],1,sum)

   # remove spaces and other characters from Sample ID
dfVirus$Sample.ID <- gsub(" ", "", dfVirus$Sample.ID, fixed = TRUE)
dfVirus$Sample.ID <- gsub("ab", "", dfVirus$Sample.ID, fixed = TRUE)
dfVirus$Sample.ID <- gsub("a+b", "", dfVirus$Sample.ID, fixed = TRUE)
dfVirus$Sample.ID <- gsub("AB", "", dfVirus$Sample.ID, fixed = TRUE)

apply(dfVirus[,allPathogens],2,sumNAs)



###############################
# Optical summary data
dfOptSum <- read.csv(file.path(raw.path,"PhaseIV","MMSDOptSummary.csv"),stringsAsFactors = FALSE,skip = 1)

# remove spaces and dashes in Sample ID
dfOptSum$FieldExpID <- gsub(" ", "", dfOptSum$FieldExpID, fixed = TRUE)
dfOptSum$FieldExpID <- gsub("-", "", dfOptSum$FieldExpID, fixed = TRUE)


dfLeftOut <- dfOptSum[!dfOptSum$FieldExpID  %in% dfmerge$Sample.ID,]

dfLeftOut2 <- dfOptSum[!dfmerge$Sample.ID %in%  dfOptSum$FieldExpID,]


# QARows <- multiGrep2(c('Blank','Replicate'),dfmerge$Comments,ignore.case = TRUE,)
# dfQA <- dfmerge[QARows,]
# df <- dfmerge[-QARows,]


########################
# Sample tracking data
dfTracking <- read.csv(file.path(raw.path,"PhaseIV","Sample tracking.csv"),stringsAsFactors = FALSE)

apply(dfTracking[,c("Virus.Sample.","Optics.DOC.Sample.", "GLWI.Sample.")],2,sum,na.rm=TRUE)

# remove spaces in Sample ID
dfTracking$Sample.ID <- gsub(" ", "", dfTracking$Sample.ID, fixed = TRUE)

# Add site abbreviations and site short names
dfTracking$abbrev <- substr(dfTracking$Sample.ID,1,2)
unique(dfTracking$abbrev)
Sites <- c('Bark','Cedar','Underwood','70th','16th')
names(Sites) <- c('BK','CG','UW','MW','MC')
dfTracking$site <- Sites[dfTracking$abbrev]
dfTracking <- dfTracking[dfTracking$abbrev %in% names(Sites),]

# Limit to samples sent to lab
samples <- apply(dfTracking[,c("Virus.Sample.","Optics.DOC.Sample.", "GLWI.Sample.")],1,sum,na.rm=TRUE)
dfTracking <- dfTracking[which(samples>0),]

###############################
# Merge all bacteria, virus, optical and sample tracking into one dataframe

#merge dfBact and dfVirus
dfmerge <- merge(dfBact, dfVirus, by="Sample.ID", all=TRUE)

dfmerge <- merge(dfmerge, dfOptSum, by.x="Sample.ID", by.y="FieldExpID", all=TRUE)

# merge mergeddf and tracking data
dfmerge <- merge(dfTracking,dfmerge,  by="Sample.ID", all=TRUE)


#################################################
#Check for consistency in data availability

# Optical
dfnoOpt <- dfmerge[is.na(dfmerge$GRnumber),]
# CG 100 and MC 100 are not reported for optical. Doesn't look like they exist.
#  -They are both on the same COC from 9/19/2012 here: M:\QW Monitoring Team\MMSD\Phase IV\virus\Optics\2012
# CG 116 and 117 had incorrect Sample IDs (GC 116 and 117). These were corrected
# MW 116 is not reported. There was no COC that had it included. Probably does not exist
# BK157 was a virus QA sample and should have indicated that there was an optical sample associated.

dfnoVirus <- dfmerge[is.na(dfmerge$ID),]
dfnoVirus2 <- dfnoVirus[which(dfnoVirus$Virus.Sample.==1),]
# One virus sample had the wrong Sample ID in the virus results. 
# It had MC 102, and it should have been MC 104. This was corrected in the 
# file with suffix "SRCRevised"
# The remaining samples here are all blank samples and are represented in a different
# results file for this project

dfnoBact <- dfmerge[is.na(dfmerge$FT),]
dfnoBact <- dfnoBact[which(dfnoBact$GLWI.Sample.==1),]
# BK 100, MW 100, and UW 100 are missing. These are blanks from the very first month of sampling
# They may not have been sent to UWM

# Check dates among the different analyses reporting
pdateTracking <- as.POSIXct(dfmerge$Start.date.time..mm.dd.yy.hh.mm.,format="%m/%d/%Y %H:%M")
pdateVirus <- as.POSIXct(dfmerge$Start.Date,format="%d-%b-%y")
pdateOpt <- as.POSIXct(paste(dfmerge$SmplDate, dfmerge$SmplStartTime,sep=" "),format="%m/%d/%Y %H:%M")
pdateBact <- as.POSIXct(paste(dfmerge$DATE, dfmerge$SmplStartTime,sep=" "),format="%m/%d/%Y %H:%M")


pdateTracking - pdateOpt #Fixed a few dates in the sample tracking that had 2-digit years: Phase IV/virus/R/Sample tracking.xls
sum(abs((pdateTracking - pdateVirus)/3600/24) > 1,na.rm=TRUE)

pdateTracking - pdateBact #Fixed a few dates in the sample tracking that had 2-digit years: Phase IV/virus/R/Sample tracking.xls
sum(abs((pdateTracking - pdateBact)/3600/24) > 7,na.rm=TRUE)
data.frame(pdateTracking,pdateBact)

############################################
# Separate QA samples

QArows <- which(dfmerge$QA.QC==1)
QArows <- c(QArows,multiGrep2(c('blank','replicate','duplicate'),dfmerge$QA,ignore.case=TRUE))
QArows <- c(QArows,multiGrep2(c('MW900','MW901','QC5','QC6','Sourcewater1'),dfmerge$Sample.ID,ignore.case=TRUE))

QArows <- unique(QArows)

dfQA <- dfmerge[QArows,]
df <- dfmerge[-QArows,]

########################################################
# Separate WW samples and Wauwatosa/Madison storm sewer samples

WWrows <- multiGrep2(c('JI','SS'),df$Sample.ID,ignore.case=TRUE)
dfWW <- df[WWrows,]
df <- df[-WWrows,]

WArows <- multiGrep2(c('WA','MAD'),df$Sample.ID,ignore.case=TRUE)
dfWA <- df[WArows,]
df <- df[-WArows,]

#########################################################################
# Check whether all virus environmental samples came through the process

dfVirusNotIncluded <- dfVirus[which(!dfVirus$Sample.ID %in% df$Sample.ID),]
sum(dfmerge$ID %in% dfVirusNotIncluded$ID)
dfCheck <- dfmerge[which(dfmerge$ID %in% dfVirusNotIncluded$ID),]
# All samples missing are QA--looks like Virus samples all made it through

#########################################################################
# Check whether all bacteria environmental samples came through the process

dfBactNotIncluded <- dfBact[which(!dfBact$Sample.ID %in% df$Sample.ID),]
sum(dfmerge$FT %in% dfBactNotIncluded$FT)
dfCheck <- dfmerge[which(dfmerge$FT %in% dfBactNotIncluded$FT),]
# All missing samples are QA. Looks like all bacteria samples made it through


#########################################################################
# Check whether all optical environmental samples came through the process

dfOptNotIncluded <- dfOptSum[which(!dfOptSum$FieldExpID %in% df$Sample.ID),]
sum(dfmerge$GRnumber %in% dfOptNotIncluded$GRnumber)
dfCheck <- dfmerge[which(dfmerge$GRnumber %in% dfOptNotIncluded$GRnumber),]
# All missing samples are QA. Looks like all optical samples made it through


# change Rdata to rds in the filenames
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write.csv(df,file=file.path(cached.path,cached.save,'VirusPhaseIVData.csv'),row.names=FALSE)
saveRDS(df,file=file.path(cached.path,cached.save,'VirusPhaseIVData.Rdata'))

write.csv(dfQA,file=file.path(cached.path,cached.save,'QAVirusPhaseIVData.csv'),row.names=FALSE)
saveRDS(dfQA,file=file.path(cached.path,cached.save,'QAVirusPhaseIVData.Rdata'))

write.csv(dfWW,file=file.path(cached.path,cached.save,'WWVirusPhaseIVData.csv'),row.names=FALSE)
saveRDS(dfWW,file=file.path(cached.path,cached.save,'WWVirusPhaseIVData.Rdata'))

write.csv(dfWA,file=file.path(cached.path,cached.save,'WAMADVirusPhaseIVData.csv'),row.names=FALSE)
saveRDS(dfWA,file=file.path(cached.path,cached.save,'WAMADVirusPhaseIVData.Rdata'))

