#Test Heuristic model parameters from MMSD with the GLRI data set

library(smwrBase)
library(caret)


plotModel <- function(m,df,response,selectedRows,colorOptions,plotColors, ...){
  predictions <- predict(m,newdata = df[selectedRows,])
  par(mfcol=c(1,1))
  plot(df[selectedRows,response],predictions,
       xlab="Observed",ylab="Predicted", col=plotColors,...)
  abline(0,1)
  mtext(paste(names(coef(m))[-1],collapse=", "),side=3,line=1,cex=0.8)
  mtext(paste("Ordinary Least Squares results:",response),side=3,line=3,font=2,cex=1)
  legend(x="topleft",legend=names(colorOptions),col=colorOptions,pch=20,text.col=colorOptions,cex=0.7)
}

plotModel_png <- function(m,df,response,selectedRows,colorOptions,plotColors, ...){
  predictions <- predict(m,newdata = df[selectedRows,])
  par(mfcol=c(1,1))
  plot(df[selectedRows,response],predictions,
       xlab="",ylab="", col=plotColors, las=1, cex=1.4,cex.axis=1.4,...)
  abline(0,1)
  mtext("Predicted (CN/100ml)",side=2,line=3,cex=1.4)
  mtext("Observed (CN/100ml)",side=1,line=3,cex=1.4)
  mtext(paste(names(coef(m))[-1],collapse=", "),side=3,line=1,cex=1.3,font=2)
  #  mtext(paste("Ordinary Least Squares results:",response),side=3,line=3,font=2,cex=1)
  legend(x="topleft",legend=names(colorOptions),col=colorOptions,pch=20,text.col=colorOptions,cex=1.2)
}

#load("GLRIdfOptSummary2016-03-03.RData") #From M:\QW Monitoring Team\Optical sensors\AqualogProcessing\R
#original bacteria data:D:\SRCSync\GLRI Toxics\Phase II\Wastewater Intensive\Data reported from labs\Merged Bacteria\GLRI03-16-16_mergedBact.RData

load("GLRIOptSummaryJan072015.RData")



df <- dfOptSumAll

df$Lachno.2.cn.100ml

df <- df[!is.na(df$Lachno.2.cn.100ml),]
df$Lachno.2.cn.100ml <- ifelse(df$Lachno.2.cn.100ml==1,112,df$Lachno.2.cn.100ml)

df <- df[!is.na(df$BACHUM.cn.100mls),]
df$BACHUM.cn.100mls <- ifelse(df$BACHUM.cn.100mls==1,112,df$BACHUM.cn.100mls)

df$loglachno <- log10(df$Lachno.2.cn.100ml)
df$logBachum <- log10(df$BACHUM.cn.100mls)
df$sinDate <- fourier(df$GMTStartTime)[,1]
df$cosDate <- fourier(df$GMTStartTime)[,2]



dummies <- dummyVars(~Site,data=df)

sites <- c("Manitowoc","Milwaukee","Menominee","Rouge","Clinton","Raisin","Portage","Maumee")
names(sites) <- c("Man","MKE","Men","RO","CL","RA","PO","MA")

df$Man <- ifelse(df$Site == "Manitowoc",1,0)
df$MKE <- ifelse(df$Site == "Milwaukee",1,0)
df$Men <- ifelse(df$Site == "Menominee",1,0)
df$RO <- ifelse(df$Site == "Rouge",1,0)
df$CL <- ifelse(df$Site == "Clinton",1,0)
df$RA <- ifelse(df$Site == "Raisin",1,0)
df$PO <- ifelse(df$Site == "Portage",1,0)
df$MA <- ifelse(df$Site == "Maumee",1,0)

df.all.sites <- df

boxplot(df.all.sites$Lachno.2.cn.100ml~df$Site,log="y",las=2)

df <- df.all.sites[df.all.sites$Site == "Manitowoc",]
df <- df.all.sites[df.all.sites$state == "WI",]
df <- df.all.sites[df.all.sites$LandUse %in% c("urban","mixed"),]
df <- df.all.sites[df.all.sites$Site %in% c("Portage","Maumee"),]
df <- df.all.sites[df.all.sites$Site %in% c("Menominee"),]
df <- df.all.sites

colorOptions <- c("orange","yellow2","skyblue","black","springgreen4","blue","grey","darkorchid1")
names(colorOptions) <- unique(df$abbrev)
selectedSiteColors <- colorOptions[unique(df$abbrev)]

m <- lm(loglachno ~  rF_T + T*sinDate + T*cosDate + F*sinDate + F*cosDate , data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

response <- "loglachno"

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)


m <- lm(logBachum ~  rF_T + T*sinDate + T*cosDate + F*sinDate + F*cosDate , data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

response <- "logBachum"

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)


# all sites with site variables
m <- lm(loglachno ~  rF_T + T*sinDate + T*cosDate + F*sinDate + F*cosDate + Man + MKE + RO + CL + RA + PO + MA, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

response <- "loglachno"

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

m <- lm(logBachum ~  rF_T + T*sinDate + T*cosDate + F*sinDate + F*cosDate + Man + MKE + RO + CL + RA + PO + MA, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

response <- "logBachum"

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

### All sites with LASSO variables plus site variables
df <- df.all.sites
m <- lm(logBachum ~   A280 + HIX_2002 + Sag275_290 + Sag350_400*sinDate + Sag350_400*cosDate + MKE + RO + CL + RA + PO + MA, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

response <- "logBachum"

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

## USE THIS IN PRESENTATION 7/24/2018 ###*****
#Stepwise from here
m2<-step(m,k=log(nrow(df)))
summary(m2)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

response <- "logBachum"

png("cached_data/2_visualize/LassoVarsAllSites.png",width=900,height=600)
plotModel_png(m=m2,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)
dev.off()

###### Most urban sites: Cl, RO
df <- df.all.sites[df.all.sites$LandUse %in% c("urban"),]

m <- lm(loglachno ~  rF_T + T*sinDate + T*cosDate + F*sinDate + F*cosDate + CL + RO, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

response <- "loglachno"

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)


m <- lm(logBachum ~  rS3.25_S1.25 + T*sinDate + T*cosDate + F*sinDate + F*cosDate + Man + MKE + RO + CL + RA + PO + MA, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

response <- "logBachum"

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

### CL and RO with LASSO variables + season + plus site variables

df <- df.all.sites[df.all.sites$LandUse %in% c("urban"),]

m <- lm(logBachum ~   A280 + HIX_2002 + Sag275_290*sinDate + Sag275_290*cosDate + Sag350_400  +  CL , data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

response <- "logBachum"

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

#Stepwise from here
m2<-step(m,k=log(nrow(df)))
summary(m2)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

response <- "logBachum"

plotModel(m=m2,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)




m <- lm(logBachum ~  Sag350_400 + A280*sinDate + A280*cosDate + Sag275_290*sinDate + Sag275_290*cosDate + CL , data = df)

m <- lm(logBachum ~  Sag350_400 + A280 + D + df$HIX_2002* + Sag275_290*sinDate + Sag275_290*cosDate + CL , data = df)

summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

response <- "logBachum"

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)


###### Most urban sites with just the low bachum values: Cl, RO
df <- df.all.sites[df.all.sites$LandUse %in% c("urban"),]
df <- subset(df,logBachum < 4)

m <- lm(logBachum ~  T*sinDate + T*cosDate + F*sinDate + F*cosDate + CL , data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

response <- "logBachum"

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)


##### OH Ag sites: 
df <- df.all.sites[df.all.sites$abbrev %in% c("MA","PO"),]
unique(df$abbrev)

m <- lm(loglachno ~  rF_T + T*sinDate + T*cosDate + F*sinDate + F*cosDate + MA + PO, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

response <- "loglachno"

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)


df <- df.all.sites[df.all.sites$abbrev %in% c("MA","PO","RM"),]
unique(df$abbrev)
selectedSiteColors <- colorOptions[unique(df$abbrev)]


m <- lm(logBachum ~  T*sinDate + T*cosDate + F*sinDate + F*cosDate + MA + PO  , data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

response <- "logBachum"

png("cached_data/2_visualize/PO_MA_T_F_season_regression.png",width=900,height=600)
plotModel_png(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)
dev.off()


##### WI Ag sites: Manitowoc has very low concentration and lots of non-detects
df <- df.all.sites[df.all.sites$abbrev %in% c("JI"),]
unique(df$abbrev)

m <- lm(loglachno ~  rF_T + T*sinDate + T*cosDate + F*sinDate + F*cosDate + MKE, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

response <- "loglachno"

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)


unique(df$abbrev)

m <- lm(logBachum ~  T*sinDate + T*cosDate + F*sinDate + F*cosDate   , data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

response <- "logBachum"

plotModel_png(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

