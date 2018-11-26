#Explore seasonality of optical signals and bacteria

library(smwrBase)
library(dataRetrieval)

#set data directories
raw.path <- "raw_data"
cached.path <- "cached_data"
summary.path <- "SummaryVariables"
summary.save <- "1_SummaryVariables"
cached.save <- "0_munge"

df <- readRDS(file.path(cached.path,summary.save,"dfOptP3P4Combined.rds"))

response <- "lachno2"
which(substr(names(df),1,1)=="A")
AbsVars <- names(df)[c(61:138,232:240)]
FlVars <- names(df)[c(17:60,139:231)]
IVs <- c(AbsVars,FlVars)

month <- as.POSIXlt(df$psdate)$mon + 1
jday <- as.POSIXlt(df$psdate)$yday + 1

boxplot(df$A~month)
boxplot(df$T~month)
boxplot(df$lachno2~month,log="y")
bp <- boxplot((df$lachno2+df$bacHum)~month,log="y")
mtext("n = ", line = 0.1, side = 3,adj=0, cex=1)
mtext(paste(bp$n, sep = ""), at = seq_along(bp$n), line = 0.1, side = 3,cex=1)

plot(log10(df$lachno2+df$bacHum)~jday)
lines(lowess(log10(df$lachno2+df$bacHum)~jday,f = 0.4))
             
fourier(df$psdate)

jdate <- as.POSIXlt(df$psdate)$yday 

jsecs <- as.numeric(df$psdate) - 
as.numeric(as.POSIXct(paste0(as.POSIXlt(df$psdate)$year+1900,"-01-01 00:00"),tz='Etc/GMT-6'))

dfQ <- readNWISuv(siteNumbers = "04087120",parameterCd = "00060",startDate = "2014-03-01",endDate = "2014-03-30")
