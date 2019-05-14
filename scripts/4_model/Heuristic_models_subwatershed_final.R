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
AbsVars <- names(df)[c(64:141, 236:246)] #define which variables are from Abs spectra 64:141, 236:246
FlVars <- names(df)[c(18:63,142:235)]   #define which variables are from Fl spectra: 17:62, 81:119
VarsSensors <- c("S1.25","S2.25","S3.25","T","C","F","A254","A295","A350","A400","CSO","UW","MC","sinDate","cosDate")
IVs <- c(grep("B",FlVars,invert = TRUE,value = TRUE))

IVs <- c(IVs,c("UW","MC","sinDate","cosDate"))

sites <- unique(df$abbrev)


selectedSites <- c("MC","MW","UW")
selectedRows <- which(df$abbrev %in% selectedSites)


colorOptions <- c("orange","yellow2","skyblue","black","springgreen4","blue","grey","darkorchid1")
names(colorOptions) <- unique(df$abbrev)
selectedSiteColors <- colorOptions[selectedSites]

df <- df[selectedRows,]
plotColors <- colorOptions[df[selectedRows,"abbrev"]]



### Interaction with stepwise for variables with only existing sensors
## Lachno2; R^2 = 0.82##
## T, F, Turbidity, sites, CSO, interaction with season

response <- "loglachno2"
#dfModel <- df[,c(response,VarsSensors[c(1,2,3,4,6,11:15)],"Turbidity_mean")]

### FINAL LACHNO MODEL  #####
dfModel <- df[,c(response,VarsSensors[c(4,6)],"UW","MC","CSO","Turbidity_mean","sinDate","cosDate")]
#dfModel <- df[,c(response,VarsSensors[c(1,3)],"UW","MC","CSO","Turbidity_mean","sinDate","cosDate")]

m_init <- lm(loglachno2 ~ ., data = dfModel)
m_step <- step(m_init, scope = . ~ .^2, direction = "both",k=log(dim(dfModel)[1]))

summary(m_step)
selectedRows <- c(1:dim(dfModel)[1])
dfModel$abbrev <- df$abbrev
plotColors <- colorOptions[dfModel[,"abbrev"]]

plotModel(m=m_step,df=dfModel,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)



### Interaction with stepwise for variables with only existing sensors
## bachuman R^2 = ##

response <- "logbacHum"
df$logbacHum <- log10(df$bacHum)
dfModel <- df[,c(response,VarsSensors[c(1,2,3,4,6,11:15)],"Turbidity_mean")]
dfModel <- df[,c(response,VarsSensors[c(4,6)],"UW","MC","CSO","Turbidity_mean","sinDate","cosDate")]

m_init <- lm(logbacHum ~ ., data = dfModel)
m_step <- step(m_init, scope = . ~ .^2, direction = "both",k=log(dim(dfModel)[1]))

summary(m_step)
selectedRows <- c(1:dim(dfModel)[1])
dfModel$abbrev <- df$abbrev
plotColors <- colorOptions[dfModel[,"abbrev"]]

plotModel(m=m_step,df=dfModel,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)


### Interaction with stepwise for variables with only existing sensors
## bachuman + Lachno R^2 = ##

response <- "logHB"
df$logHB <- log10(df$bacHum + df$lachno2)
dfModel <- df[,c(response,VarsSensors[c(1,2,3,4,6,11:15)],"Turbidity_mean")]
dfModel <- df[,c(response,VarsSensors[c(4,6)],"UW","MC","CSO","Turbidity_mean","sinDate","cosDate")]

m_init <- lm(logHB ~ ., data = dfModel)
m_step <- step(m_init, scope = . ~ .^2, direction = "both",k=log(dim(dfModel)[1]))

summary(m_step)
selectedRows <- c(1:dim(dfModel)[1])
dfModel$abbrev <- df$abbrev
plotColors <- colorOptions[dfModel[,"abbrev"]]

plotModel(m=m_step,df=dfModel,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

## Final HB model ##
form <- formula(logHB ~ T + F*sinDate + F*cosDate + Turbidity_mean*sinDate + Turbidity_mean*cosDate + UW + MC + CSO)
m_HB <- lm(form, data = df)
summary(m_HB)

##  ##  ##  ##


### Interaction with stepwise for variables with only existing sensors
## E coli R^2 = ##

response <- "logEcoli"
dfModel <- df[,c(response,VarsSensors[c(1,2,3,4,6,11:15)],"Turbidity_mean","Water_Temperature_mean")]
dfModel <- df[,c(response,VarsSensors[c(4,6)],"UW","MC","CSO","Turbidity_mean","sinDate","cosDate")]

m_init <- lm(logEcoli ~ ., data = dfModel)
m_step <- step(m_init, scope = . ~ .^2, direction = "both",k=log(dim(dfModel)[1]))

summary(m_step)
selectedRows <- c(1:dim(dfModel)[1])
dfModel$abbrev <- df$abbrev
plotColors <- colorOptions[dfModel[,"abbrev"]]

plotModel(m=m_step,df=dfModel,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

## Final  E coli model ##
form <- formula(logEcoli ~ T + F*sinDate + F*cosDate + Turbidity_mean*sinDate + Turbidity_mean*cosDate + UW + MC + CSO)
m_Ecoli <- lm(form, data = df)
summary(m_Ecoli)

##  ##  ##  ##


### Interaction with stepwise for variables with only existing sensors
## Enterococci R^2 = ##

response <- "logEnt"
df$logEnt <- log10(df$ent)
dfModel <- df[,c(response,VarsSensors[c(1,2,3,4,6,11:15)],"Turbidity_mean","Water_Temperature_mean")]
dfModel <- df[,c(response,VarsSensors[c(4,6)],"UW","MC","CSO","Turbidity_mean","sinDate","cosDate")]

m_init <- lm(logEnt ~ ., data = dfModel)
m_step <- step(m_init, scope = . ~ .^2, direction = "both",k=log(dim(dfModel)[1]))

summary(m_step)
selectedRows <- c(1:dim(dfModel)[1])
dfModel$abbrev <- df$abbrev
plotColors <- colorOptions[dfModel[,"abbrev"]]

plotModel(m=m_step,df=dfModel,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

## Final  E coli model ##
form <- formula(logEnt ~ T*sinDate + T*cosDate + F*sinDate + F*cosDate + Turbidity_mean*sinDate + Turbidity_mean*cosDate + UW + MC + CSO)
m_Ent <- lm(form, data = df)
summary(m_Ent)

##  ##  ##  ##