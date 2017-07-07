#LassoOpticalRegressions

library(glmnet)
library(smwrBase)
library(survival)
library(MASS)

#set data directories
raw.path <- "raw_data"
cached.path <- "cached_data"
summary.path <- "SummaryVariables"
summary.save <- "1_SummaryVariables"
cached.save <- "0_munge"

df <- readRDS(file.path(cached.path,summary.save,"dfOptP3P4Combined.rds"))

response <- "eColi"
which(substr(names(df),1,1)=="A")
AbsVars <- names(df)[c(61:138,232:240)]
FlVars <- names(df)[c(17:60,139:231)]
write.csv(FlVars,file=file.path(cached.path,summary.save,"FlVars.csv"),row.names = FALSE)
dfSensors <- read.csv(file=file.path(cached.path,summary.save,"FlVarsSensors.csv"), stringsAsFactors = FALSE)
sensorVars <- dfSensors[which(dfSensors$Sensor==1),"Variable"]
IVs <- c(sensorVars,c("UW","MC","sinDate","cosDate"))
IVs <- c(sensorVars,c("sinDate","cosDate"))
penalty.factor <- c(rep(1,(length(IVs)-2)),0,0)

sites <- unique(df$abbrev)

selectedRows <- which(df$abbrev %in% c("HW","MC","MW","UW","MF"))
selectedRows <- which(df$abbrev %in% c("MC","MW","UW"))
selectedRows <- which(df$abbrev %in% sites)
selectedRows <- which(df$abbrev %in% c("HW","MC","MW","UW","MF","CG","LD"))
selectedRows <- which(df$abbrev %in% c("HW","MC","MW","UW","MF","CG","LD","BK"))
selectedRows <- which(df$abbrev %in% c("MC","UW","MW"))

x <- as.matrix(df[selectedRows,IVs])
y <- log10(df[selectedRows,response])

glm.family <- "gaussian"
#glm.family <- "poisson"

foldid <- as.numeric(as.factor(df[selectedRows,"abbrev"])) #use sites as folds
mg.cv <- cv.glmnet(x=x, y=y,foldid=foldid,family=glm.family)
#mg.cv <- cv.glmnet(x=x, y=y,nfolds = 10, family=glm.family, penalty.factor=penalty.factor)
mg <- glmnet(x=x, y=y,family=glm.family, penalty.factor=penalty.factor)

# #use lambda with min cross-validated error
# mg.minl <- glmnet(x=x, y=y, s = mg.cv$lambda.min,family=glm.family)

#Extract Coefficients from cv-determined model
Coefficients <- coef(mg, s = mg.cv$lambda.1se)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index];Active.Coefficients
Active.Coef.names <- row.names(coef(mg.cv))[Active.Index];Active.Coef.names

#Plot cross validated errors and other model results
plot(mg.cv)
summary(mg.cv)


#Plot Observed vs Predicted from cv-determined model

colorOptions <- c("orange","yellow2","skyblue","black","springgreen4","blue","grey","darkorchid1")
names(colorOptions) <- unique(df$abbrev)

plotColors <- colorOptions[df[selectedRows,"abbrev"]]
levels(as.factor(df[selectedRows,"abbrev"]))

axis.log <- ""
axis.limits <- c(1,7)

par(mfcol=c(1,1))
plot(y,predict(mg.cv,newx=x,s="lambda.1se",type="response"),
     xlab="Observed",ylab="Predicted",col=plotColors,pch=20,log=axis.log,ylim=axis.limits,xlim=axis.limits)
abline(0,1)
mtext(paste(Active.Coef.names[2:length(Active.Coef.names)],collapse=", "),side=3,line=1,cex=0.8)
mtext(paste("Cross-validated Lasso model for ",response),side=3,line=3,font=2,cex=1)
mtext(paste(glm.family," GLM regression"),side=3,line=2,font=2,cex=1)
legend(x="topleft",legend=names(colorOptions),col=colorOptions,pch=20,text.col=colorOptions,cex=0.7)


#Use Lasso-determined variables in tobit regression for censored values. Use
#stepwise regression to narrow down the variables further.
Active.Coef.names
if("sinDate" %in% Active.Coef.names & !("cosDate" %in% Active.Coef.names)){
   Active.Coef.names <- c(Active.Coef.names,"cosDate")}
if("cosDate" %in% Active.Coef.names & !("sinDate" %in% Active.Coef.names)){
  Active.Coef.names <- c(Active.Coef.names,"sinDate")}

#Tobit regression
y.orig <- log10(df[selectedRows,response])
y <- Surv(y.orig, df[selectedRows,response]>225, type="left")
dfPredStd <- as.data.frame(df[selectedRows,IVs])
formUpper <- formula(paste('y ~',paste(Active.Coef.names[-1],collapse=' + ')))
formLower <- formula(paste('y ~',paste(c("sinDate","cosDate"),collapse=' + ')))
msurvStd <- survreg(formUpper,data=dfPredStd,dist='weibull')
summary(msurvStd)

#Stepwise regresson using tobit results as starting point
msurvStep <- stepAIC(msurvStd,list(upper=formUpper,lower=formLower))
summary(msurvStep)

#Plot stepwise results
predictions <- predict(msurvStep,newdata = df[selectedRows,IVs])
par(mfcol=c(1,1))
plot(y.orig,predictions,
     xlab="Observed",ylab="Predicted",col=plotColors,pch=20,log=axis.log,ylim=axis.limits,xlim=axis.limits)
abline(0,1)
mtext(paste(names(coef(msurvStep))[-1],collapse=", "),side=3,line=1,cex=0.8)
mtext(paste("Tobit from Lasso model for ",response),side=3,line=3,font=2,cex=1)
mtext(paste(glm.family," GLM regression"),side=3,line=2,font=2,cex=1)
legend(x="topleft",legend=names(colorOptions),col=colorOptions,pch=20,text.col=colorOptions,cex=0.7)



plotModel <- function(m,df,response,selectedRows,colorOptions,plotColors, ...){
  #Plot stepwise results
  predictions <- predict(m,newdata = df[selectedRows,])
  par(mfcol=c(1,1))
  plot(df[selectedRows,response],predictions,
       xlab="Observed",ylab="Predicted", col=plotColors,...)
  abline(0,1)
  mtext(paste(names(coef(m))[-1],collapse=", "),side=3,line=1,cex=0.8)
  mtext(paste(i,"-variable model from Tobit/Lasso results",response),side=3,line=3,font=2,cex=1)
  legend(x="topleft",legend=names(colorOptions),col=colorOptions,pch=20,text.col=colorOptions,cex=0.7)
}

variables <- names(coef(msurvStep))[-1]


selectedRows <- which(!is.na(df[,response]))
filenm <- "plotStepsSensorVars3SitesEColi.pdf"
pdf(filenm)
for(i in 1:length(variables)){
form <- formula(paste('y ~',paste(c(variables[1:i]),collapse=' + ')))
y <- Surv(log10(df[selectedRows,response]), df[selectedRows,response]>225, type="left")
m <- survreg(form,data=df,dist='weibull')
plotColors <- colorOptions[df[selectedRows,"abbrev"]]
plotModel(m=m,df=df,response=response,selectedRows=selectedRows,
          colorOptions=colorOptions,plotColors=plotColors,
          pch=20,log="x")
}
dev.off()
shell.exec(filenm)


#For the Gaussian regression, transform back to real space and graph with log axis
#to be comparable. Use smearing coefficient to correct for bias

# mg.residuals <- predict(mg,newx=IVmatrix, s = mg.cv$lambda.min,type="response")-df[,responsevar]
# smear.coef <- sum(10^(mg.residuals))/length(mg.residuals)
# LD <- 10^(predict(mg,newx=IVmatrix, s = mg.cv$lambda.min,type="response"))*smear.coef
# 

