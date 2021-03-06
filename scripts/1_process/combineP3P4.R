# Merge summary optical data files from Phase III and Phase IV efforts

library(smwrBase) # Use this to determine seasonality terms

#set data directories
raw.path <- "raw_data"
cached.path <- "cached_data"
summary.path <- "SummaryVariables"
summary.save <- "1_SummaryVariables"
cached.save <- "0_munge"

dfP3 <- readRDS(file.path(cached.path,summary.save,"dfOptSummaryMMSDP3.rds"))
dfP4 <- readRDS(file.path(cached.path,summary.save,"dfOptSummaryMMSDP4.rds"))

dfP4$psdate <- as.POSIXct(dfP4$Start.date.time..mm.dd.yy.hh.mm., format= "%m/%d/%Y %H:%M", tz = 'Etc/GMT+6')
dfP4$pedate <- as.POSIXct(dfP4$End.date.time..mm.dd.yy.hh.mm., format= "%m/%d/%Y %H:%M", tz = 'Etc/GMT+6')
dfP3$psdate <- as.POSIXct(format(as.POSIXct(dfP3$psdate),tz='Etc/GMT+6',usetz=TRUE),tz='Etc/GMT+6')
dfP3$pedate <- as.POSIXct(format(as.POSIXct(dfP3$pedate),tz='Etc/GMT+6',usetz=TRUE),tz='Etc/GMT+6')

#names(dfP3)[which(names(dfP3) == "Abb")]

#Change names to coincide between P3 and P4
names(dfP3)[names(dfP3) %in% "Abb"] <- "abbrev"
names(dfP3)[names(dfP3) %in% "BacHumancen"] <- "bacHum"
names(dfP3)[names(dfP3) %in% "lachnocen"] <- "lachno2"
names(dfP3)[names(dfP3) %in% "Entero.CN.100ml"] <- "ent"
names(dfP3)[names(dfP3) %in% "E..coli.CFU.100ml"] <- "eColi"
names(dfP3)[names(dfP3) %in% "Human_virus"] <- "humanVirus"
names(dfP3)[names(dfP3) %in% "GLWI.FT"] <- "FT"

dfP3$abbrev <- as.character(dfP3$abbrev)
dfP3$OMabs_dilution <- as.character(dfP3$OMabs_dilution)

dfP3$hydro_condition <- as.character(dfP3$Sample_Event)
dfP4$hydro_condition <- Hmisc::capitalize(dfP4$hydro_condition)
names(dfP4)[names(dfP4) %in% "lachno"] <- "lachno2"


commonNames <- names(dfP3)[which(names(dfP3) %in% names(dfP4))]

begin <- 1
end <- length(commonNames)
df <- rbind(dfP3[,commonNames[begin:end]],dfP4[,commonNames[begin:end]])

#Add seasonal and site terms
df$sinDate <- fourier(df$psdate)[,1]
df$cosDate <- fourier(df$psdate)[,2]
df$BK <- ifelse(df$abbrev=="BK",1,0)
df$MF <- ifelse(df$abbrev=="MF",1,0)
df$UW <- ifelse(df$abbrev=="UW",1,0)
df$LD <- ifelse(df$abbrev=="LD",1,0)
df$HW <- ifelse(df$abbrev=="HW",1,0)
df$MW <- ifelse(df$abbrev=="MW",1,0)
df$MC <- ifelse(df$abbrev=="MC",1,0)
df$CG <- ifelse(df$abbrev=="CG",1,0)

c("BK","MF","UW","LD","HW","MW","MC","CG")

saveRDS(df,file=file.path(cached.path,summary.save,"dfOptP3P4Combined.rds"))

# write.csv(names(dfP3),file = "P3names.csv",row.names=FALSE)
# write.csv(names(dfP4),file = "P4names.csv",row.names=FALSE)
