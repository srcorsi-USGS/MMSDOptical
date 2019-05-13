#Test regressions with common optical signals

library(glmnet)
library(smwrBase)
library(survival)
library(MASS)

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
       xlab="",ylab="", col=plotColors, las=1, cex.axis=1.4,...)
  abline(0,1)
  mtext("Predicted (CN/100ml)",side=2,line=3,cex=1.4)
  mtext("Observed (CN/100ml)",side=1,line=3,cex=1.4)
  mtext(paste(names(coef(m))[-1],collapse=", "),side=3,line=1,cex=1.3,font=2)
#  mtext(paste("Ordinary Least Squares results:",response),side=3,line=3,font=2,cex=1)
  legend(x="topleft",legend=names(colorOptions),col=colorOptions,pch=20,text.col=colorOptions,cex=1.2)
}

#set data directories
raw.path <- "raw_data"
cached.path <- "cached_data"
summary.path <- "SummaryVariables"
summary.save <- "1_SummaryVariables"
cached.save <- "0_munge"

df <- readRDS(file.path(cached.path,summary.save,"dfOptP3P4CombinedContWQ.rds"))
df.orig <- df
DA <- c(89.9, 20.7, 26.7, 46.9, 319, 378)
names(DA) <- c("MF", "LD", "HW", "UW", "MW", "MC")
df$DA <- DA[df$abbrev]
df$Discharge_mean_DA <- df$Discharge_mean/df$DA
df$Discharge_max_DA <- df$Discharge_max/df$DA

df$CSOVol <- 0
df[which(df$GRnumber=="gr15841"),"CSOVol"] <- 594.7
df[which(df$GRnumber=="gr15887"),"CSOVol"] <- 524.9
df[which(df$GRnumber=="gr17417"),"CSOVol"] <- 341.2
df$CSO <- ifelse(df$CSOVol > 0,1,0)
df$loglachno2 <- log10(df$lachno2)
df$logEcoli <- log10(df$eColi)

response <- "lachno2"
df$logResponse <- log10(df[,response])
response <- paste0("log",response)
names(df)[dim(df)[2]] <- response

which(substr(names(df),1,1)=="A")
AbsVars <- names(df)[c(61:138,232:240)] #define which variables are from Abs spectra 63:80, 120:130
FlVars <- names(df)[c(17:60,139:231)]   #define which variables are from Fl spectra: 17:62, 81:119
IVs <- c(grep("B",FlVars,invert = TRUE,value = TRUE))

IVs <- c(IVs,c("UW","MC","sinDate","cosDate"))
#penalty.factor <- c(rep(1,length(c(AbsVars,FlVars))),0,0,0,0)

sites <- unique(df$abbrev)


#selectedSites <- c("HW","MC","MW","UW")
selectedSites <- c("MC","MW","UW")
#selectedSites <- sites
#selectedSites <- c("HW","MC","MW","UW","MF","CG","LD")
#selectedSites <- c("HW","MC","MW","UW","MF","LD")
selectedRows <- which(df$abbrev %in% selectedSites)

x <- as.matrix(df[selectedRows,IVs])
y <- log10(df[selectedRows,response])
#Plot Observed vs Predicted from cv-determined model

colorOptions <- c("orange","yellow2","skyblue","black","springgreen4","blue","grey","darkorchid1")
names(colorOptions) <- unique(df$abbrev)
selectedSiteColors <- colorOptions[selectedSites]

df <- df[selectedRows,]
plotColors <- colorOptions[df[selectedRows,"abbrev"]]


m <- lm(log10(lachno2)~A440*F1 + A440*sinDate + A440*cosDate + F1*sinDate + F1*cosDate, data = df)
summary(m)
pred <- predict(m,newdata = df)
plotColors <- colorOptions[df[,"abbrev"]]
plot(df$lachno2~pred,log="y",col=plotColors,pch=20)

m <- lm(log10(eColi)~A440*F + A440*sinDate + A440*cosDate + F*sinDate + F*cosDate, data = df)
summary(m)
pred <- predict(m,newdata = df)
plotColors <- colorOptions[df[,"abbrev"]]
plot(df$eColi~pred,log="y",col=plotColors,pch=20)

m <- lm(log10(lachno2)~T*sinDate + T*cosDate, data = df)
summary(m)
pred <- predict(m,newdata = df)
plotColors <- colorOptions[df[,"abbrev"]]
plot(df$lachno2~pred,log="y",col=plotColors,pch=20)

m <- lm(log10(bacHum)~T*sinDate + T*cosDate, data = df)
summary(m)
pred <- predict(m,newdata = df)
plotColors <- colorOptions[df[,"abbrev"]]
plot(df$bacHum~pred,log="y",col=plotColors,pch=20)

m <- lm(log10(lachno2)~rF_T + T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW, data = df)
summary(m)
pred <- predict(m,newdata = df)
plotColors <- colorOptions[df[,"abbrev"]]
plot(df$lachno2~pred,log="y",col=plotColors,pch=20)
legend(x="topleft",legend=names(selectedSiteColors),col=selectedSiteColors,pch=20,text.col=colorOptions,cex=0.7)
abline(0,1,col="red",lty=2)

m <- lm(log10(lachno2)~ CSO + Turbidity_mean + Water_Temperature_mean*Discharge_mean_DA + Water_Temperature_mean*Discharge_mean_DA +  MC + UW, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

m <- lm(log10(lachno2)~ CSO + Turbidity_mean + sinDate*Discharge_mean_DA + cosDate*Discharge_mean_DA +  MC + UW, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)




# 1-a. Potential final model, but it contains the CSO variable

filenm <- "MMSDOpticalRegression_3Sites20180723.png"
png(filenm,width = 900,height = 600)
m <- lm(log10(lachno2)~CSO + T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel_png(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20,cex=1.3)
dev.off()
shell.exec(filenm)

#Alternative with flow at 16th interation with date **** TRY USING AN ANTECEDENT BASEFLOW VARIABLE FOR 16TH.
m <- lm(log10(lachno2)~ (Discharge_max_DA*MC):sinDate + (Discharge_max_DA*MC):cosDate + T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW, data = df)
summary(m)

selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)


m <- lm(log10(lachno2)~ rF_T*sinDate + rF_T*cosDate + T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)


# 1-b. Potential final model without the CSO data
dfb <- subset(df,lachno2 < 200000)

m <- lm(log10(lachno2)~ T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW, data = dfb)
summary(m)
selectedRows <- c(1:dim(dbf)[1])
plotColors <- colorOptions[dfb[,"abbrev"]]

plotModel(m=m,df=dfb,response=response,selectedRows = c(1:dim(dfb)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

m <- lm(log10(lachno2)~ T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW, data = df)

selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

m <- lm(log10(lachno2)~ rF_T*sinDate + rF_T*cosDate + T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

# 2. Potential final model, but it contains the CSO variable

m <- lm(log10(lachno2)~CSO + Turbidity_mean + rS2.25_S3.25*sinDate + rS2.25_S3.25*cosDate + S3.25*sinDate + S3.25*cosDate + MC + UW, data = df)
summary(m)

selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

# 3. without the CSO variable

m <- lm(log10(lachno2)~ Turbidity_mean + rS2.25_S3.25*sinDate + rS2.25_S3.25*cosDate + S3.25*sinDate + S3.25*cosDate + MC + UW, data = df)
summary(m)

selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

# 4. Add water quality variables to model #1

m <- lm(log10(lachno2)~ Turbidity_mean  + CSO + T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)



m <- lm(log10(lachno2)~ Turbidity_mean  + CSO + T*Water_Temperature_mean + F*Water_Temperature_mean + MC + UW, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)


###############

m <- lm(log10(lachno2)~ Turbidity_mean  + CSO + S1.25*sinDate + S1.25*cosDate + S3.25*sinDate + S3.25*cosDate + MC + UW, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

### Test water quality without optical
m <- lm(log10(lachno2)~  CSO + Turbidity_mean*sinDate + Turbidity_mean*cosDate + MC + UW, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)


m <- lm(log10(lachno2)~ CSO + Turbidity_mean + T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)


m <- lm(log10(lachno2)~ A440*sinDate + A440*cosDate + T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW, data = df)
summary(m)

selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)


test <- df[which(df$lachno2 > 3e+05),]


m <- lm(log10(lachno2)~ df$Sag290_350 + df$B + df$HIX_2002*sinDate + T*cosDate, data = df)
summary(m)
pred <- predict(m,newdata = df)
plot(df$bacHum~pred,log="y",col=plotColors,pch=20)


# E Coli model

df$logEcoli <- log10(df$eColi)
response <- "logEcoli"

m <- lm(logEcoli ~ T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

