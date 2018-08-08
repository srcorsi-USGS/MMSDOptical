#Test Heuristic model parameters from MMSD with the GLRI data set

library(smwrBase)
library(caret)
library(glmnet)


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

dfIVs <- read.csv("GLRIIVs.csv",stringsAsFactors = FALSE)


df$Lachno.2.cn.100ml

df <- df[!is.na(df$Lachno.2.cn.100ml),]
df$Lachno.2.cn.100ml <- ifelse(df$Lachno.2.cn.100ml==1,112,df$Lachno.2.cn.100ml)

df <- df[!is.na(df$BACHUM.cn.100mls),]
df$BACHUM.cn.100mls <- ifelse(df$BACHUM.cn.100mls==1,112,df$BACHUM.cn.100mls)

df$loglachno <- log10(df$Lachno.2.cn.100ml)
df$logBachum <- log10(df$BACHUM.cn.100mls)
df$sinDate <- fourier(df$GMTStartTime)[,1]
df$cosDate <- fourier(df$GMTStartTime)[,2]

IVvars <- dfIVs$IV[which(dfIVs$Keep==1)]
IVs <- c(paste0(IVvars,"*sinDate"),paste0(IVvars,"*cosDate"))

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

dfAll <- df

boxplot(df.all.sites$Lachno.2.cn.100ml~df$Site,log="y",las=2)

#df <- df.all.sites[df.all.sites$Site == "Manitowoc",]
#df <- df.all.sites[df.all.sites$state == "WI",]
dfAll <- df[df$LandUse %in% c("urban"),]
unique(df$Site)

response <- "logBachum"
i <- 1

df <- dfAll[,c("Site",response[i],IVvars,"CL","sinDate","cosDate")]
df <- na.omit(df)
siteColors <- c("blue","green")
names(siteColors) <- c("Clinton","Rouge")

df$plotColors <- siteColors[df[,"Site"]]
IVmatrix <- as.matrix(df[,IVvars])
df$response <- (df[,response[i]])
responsevar <- "response"

##############  LASSO ###################
glm.family="gaussian"
form <- formula(paste("response ~ ",paste(c(IVs,"CL"),collapse = "+"))  )
dataMatrix <- model.matrix(form,data=df)


mg.cv <- cv.glmnet(x=IVmatrix, y=df[,responsevar],family=glm.family,nfolds=5)#,offset=fitted(msite
#mg.cv <- cv.glmnet(x=IVmatrix, y=df[,responsevar],family=glm.family,nfolds=5)
mg <- glmnet(x=IVmatrix, y=df[,responsevar],family=glm.family)


#use lambda with min cross-validated error
#mg.minl <- glmnet(x=IVmatrix, y=df[,responsevar], s = mg.cv$lambda.min,family=glm.family)

#Extract Coefficients from cv-determined model
Coefficients <- coef(mg, s = mg.cv$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index];Active.Coefficients
Active.Coef.names <- row.names(coef(mg.cv))[Active.Index];Active.Coef.names
observed <- df[,responsevar]
predicted <- as.numeric(predict(mg,newx=IVmatrix, s = mg.cv$lambda.min,type="response"))
#Plot cross validated errors and other model results
plot(mg.cv)
mtext(paste("LASSO Xval for",response[i]))
plot(observed,predicted)

png("cached_data/2_visualize/urbanLASSOplot.png",width=900,height=600)
par(mfcol=c(1,1))
plot(observed,predicted,
     xlab="",ylab="",col=df$plotColors,pch=20,log="",cex=1.4,las=1)
abline(0,1)
mtext(paste(Active.Coef.names[2:length(Active.Coef.names)],collapse=", "),side=3,line=1,cex=0.8)
mtext(paste("Cross-validated Lasso model for ",response),side=3,line=3,font=2,cex=1)
mtext(paste(glm.family," GLM regression"),side=3,line=2,font=2,cex=1)
mtext("Observed",side=1,line=3,cex=1.4)
mtext("Predicted",side=2,line=3,cex=1.4)
legend(x="topleft",legend=c("Clinton","Rouge"),col=c("blue","green"),pch=20,text.col=c("blue","green"),cex=1.4)
dev.off()

