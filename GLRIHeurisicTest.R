#Test Heuristic model parameters from MMSD with the GLRI data set

library(smwrBase)
library(caret)
load("GLRIdfOptSummary2016-03-03.RData") #From M:\QW Monitoring Team\Optical sensors\AqualogProcessing\R
                                         #original bacteria data:D:\SRCSync\GLRI Toxics\Phase II\Wastewater Intensive\Data reported from labs\Merged Bacteria\GLRI03-16-16_mergedBact.RData

#load("GLRIOptSummaryJan072015.RData")


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



df <- dfOptSumAll

df$Lachno.2.cn.100ml

df <- df[!is.na(df$Lachno.2.cn.100ml),]

df$loglachno <- log10(df$Lachno.2.cn.100ml)
df$sinDate <- fourier(df$GMTStartTime)[,1]
df$cosDate <- fourier(df$GMTStartTime)[,2]



dummies <- dummyVars(~Site,data=df)

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

m <- lm(loglachno ~  rF_T + T*sinDate + T*cosDate + F*sinDate + F*cosDate , data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

response <- "loglachno"

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)
