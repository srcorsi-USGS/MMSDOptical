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

#set data directories
raw.path <- "raw_data"
cached.path <- "cached_data"
summary.path <- "SummaryVariables"
summary.save <- "1_SummaryVariables"
cached.save <- "0_munge"

df <- readRDS(file.path(cached.path,summary.save,"dfOptP3P4CombinedContWQ.rds"))
df.orig <- df
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
AbsVars <- names(df)[c(61:138,232:240)]
FlVars <- names(df)[c(17:60,139:231)]
IVs <- c(grep("B",FlVars,invert = TRUE,value = TRUE))

IVs <- c(IVs,c("UW","MC","sinDate","cosDate"))
penalty.factor <- c(rep(1,length(c(AbsVars,FlVars))),0,0,0,0)

sites <- unique(df$abbrev)

selectedSites <- c("HW","MC","MW","UW")
selectedSites <- c("MC","MW","UW")
selectedSites <- sites
selectedSites <- c("HW","MC","MW","UW","MF","CG","LD")
selectedSites <- c("HW","MC","MW","UW","MF","LD")
selectedSites <- c("MC","UW","MW")
selectedRows <- which(df$abbrev %in% selectedSites)

x <- as.matrix(df[selectedRows,IVs])
y <- log10(df[selectedRows,response])
#Plot Observed vs Predicted from cv-determined model

colorOptions <- c("orange","yellow2","skyblue","black","springgreen4","blue","grey","darkorchid1")
names(colorOptions) <- unique(df$abbrev)
selectedSiteColors <- colorOptions[selectedSites]

df <- df[selectedRows,]
plotColors <- colorOptions[df[selectedRows,"abbrev"]]


# 4. Add water quality variables to model #1

m <- lm(log10(lachno2)~ Turbidity_mean  + CSO + PT*sinDate + PT*cosDate + PA*sinDate + PA*cosDate + MC + UW, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)



m <- lm(log10(lachno2)~ Turbidity_mean  + CSO + PT*Water_Temperature_mean + PA*Water_Temperature_mean + MC + UW, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)




m <- lm(log10(lachno2)~ Turbidity_mean  + CSO + PT*Water_Temperature_mean + PA*Water_Temperature_mean + MC + UW, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)


m <- lm(log10(lachno2)~ Turbidity_mean  + CSO + PT*Water_Temperature_mean + F*Water_Temperature_mean + MC + UW, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)


