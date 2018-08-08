library(glmnet)
library(leaps)

load("dfAll.RData")
siteColors <- colors()[c(24,498,589,614,6,26,585,96)]
names(siteColors) <- unique(dfAll$Site)

dfClasses <- read.csv("dfAllNames.csv",as.is=TRUE)
dfClasses <- read.csv(paste(Project,"/dfAllNames.csv",sep=""),as.is=TRUE)
IVvars <- dfClasses$variables[dfClasses$group=="Optical"]
IVvars <- IVvars[which(IVvars=="TDNResult"):length(IVvars)]

selectSites <- names(siteColors)
#selectSites <- c("Milwaukee","Rouge","Clinton")
#selectSites <- c("Clinton")
siteColors <- siteColors[selectSites]
responsevar <- "Lachno.2.cn.100ml"
df <- dfAll[,c("Site",responsevar,IVvars)]
df <- df[df$Site %in% selectSites,]
df <- na.omit(df)
df$plotColors <- siteColors[df[,"Site"]]
IVmatrix <- as.matrix(df[,IVvars])
df$response <- log10(df[,responsevar])
responsevar <- "response"

##############  LASSO ###################
glm.family="gaussian"
mg.cv <- cv.glmnet(x=IVmatrix, y=df[,responsevar],family=glm.family,nfolds=5)#,offset=fitted(msite),foldid=df$foldid)
mg <- glmnet(x=IVmatrix, y=df[,responsevar],family=glm.family)

#use lambda with min cross-validated error
mg.minl <- glmnet(x=IVmatrix, y=df[,responsevar], s = mg.cv$lambda.min,family=glm.family)

#Extract Coefficients from cv-determined model
Coefficients <- coef(mg, s = mg.cv$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index];Active.Coefficients
Active.Coef.names <- row.names(coef(mg.cv))[Active.Index];Active.Coef.names
observed <- df[,responsevar]
predicted <- as.numeric(predict(mg,newx=IVmatrix, s = mg.cv$lambda.min,type="response"))
#Plot cross validated errors and other model results
plot(mg.cv)
plot(observed,predicted,col=siteColors,pch=20)
abline(0,1,lty=2)
###################################################
## Use lasso results to run regsubsets with BIC
dfmodel <- data.frame(response=df[,responsevar])
names(dfmodel) <- responsevar

form <- formula(paste(responsevar,"~",paste(Active.Coef.names[-1],collapse="+")))
mregsub <- regsubsets(x=form,data=df)
summary(mregsub)
summary(mregsub)$which
summary(mregsub)$adjr2
summary(mregsub)$rsq
summary(mregsub)$cp
summary(mregsub)$obj
summary(mregsub)$bic

mregsub
regOpt <- which.min(summary(mregsub)$bic)
regOpt <- 4
regOptActive <- names(coef(mregsub,regOpt)[-1])
form <- formula(paste(responsevar,"~",paste(regOptActive,collapse="+")))
# form <- formula("responsevar ~ P_Mex300em390 + P_Sag275and295 + P_Sag350and400 + P_HIOhno2008")
# form <- formula("responsevar ~ P_Mex300em390 + P_Sag275and295 + P_Sag350and400")
# form <- formula("responsevar ~ P_fDOMex370em460  + P_Sag275and295 + P_HIOhno2008")

# form <- formula("responsevar ~ P_fDOMex370em460  + P_Sag350and400 + P_HIOhno2008")
# form <- formula("responsevar ~ P_fDOMex370em460")

mOpt <- lm(formula=form,data=df)
lmActive <- names(summary(mOpt)$coef[,1])[-1]
adjR2Opt <- summary(mOpt)$adj.r.squared

predicted <- predict(mOpt)
##
##################################################

axis.limits<- range(c(observed,predicted))

png("LassoLachnoUrban.png",,width=1440,height=960,pointsize=30)
plot(observed,predicted,ylim=axis.limits,xlim=axis.limits,col=df$plotColors,
     main="",pch=20,xlab="",ylab="")
abline(0,1,lty=2)
legend("topleft",legend=names(siteColors),col=siteColors,pch=20,cex=0.8)
mtext("Lachnospiraceae 2 Preliminary Regression with Optical Properties",side=3,line=3,font=2)
mtext(paste("adjR^2 =",round(adjR2Opt,2)),side=3,line=2,font=1)
mtext("Log10(Lachno 2 predicted)",side=2,line=3,font=2)
mtext("Log10(Lachno 2 observed)",side=1,line=3,font=2)
parmnames <- sub("P_","",lmActive)
mtext(paste(parmnames,collapse="; "),line=0.5)
dev.off()
shell.exec("LassoLachnoUrban.png")


#################################################################################
#################################################################################
# LASSO for Bacteroides
siteColors <- colors()[c(24,498,589,614,6,26,585,96)]
names(siteColors) <- unique(dfAll$Site)
IVvars <- dfClasses$variables[dfClasses$group=="Optical"]
IVvars <- IVvars[which(IVvars=="A254"):length(IVvars)]

selectSites <- c("Milwaukee","Rouge","Clinton")
#selectSites <- c("Clinton")

siteColors <- siteColors[selectSites]
#dfAll$E..coli.CFUs.100ml
responsevar <- "E..coli.CFUs.100ml"
df <- dfAll[,c("Site",responsevar,IVvars)]
df <- df[df$Site %in% selectSites,]

df <- na.omit(df)
df$plotColors <- siteColors[df[,"Site"]]
IVmatrix <- as.matrix(df[,IVvars])
df$response <- log10(df[,responsevar])
responsevar <- "response"

##############  LASSO ###################
glm.family="gaussian"
mg.cv <- cv.glmnet(x=IVmatrix, y=df[,responsevar],family=glm.family,nfolds=5)#,offset=fitted(msite),foldid=df$foldid)
mg <- glmnet(x=IVmatrix, y=df[,responsevar],family=glm.family)

#use lambda with min cross-validated error
mg.minl <- glmnet(x=IVmatrix, y=df[,responsevar], s = mg.cv$lambda.min,family=glm.family)

#Extract Coefficients from cv-determined model
Coefficients <- coef(mg, s = mg.cv$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index];Active.Coefficients
Active.Coef.names <- row.names(coef(mg.cv))[Active.Index];Active.Coef.names
observed <- df[,responsevar]
predicted <- as.numeric(predict(mg,newx=IVmatrix, s = mg.cv$lambda.min,type="response"))
#Plot cross validated errors and other model results
#plot(mg.cv)
#plot(observed,predicted)

###################################################
## Use lasso results to run regsubsets with BIC
dfmodel <- data.frame(response=df[,responsevar])
names(dfmodel) <- responsevar

form <- formula(paste(responsevar,"~",paste(Active.Coef.names[-1],collapse="+")))
mregsub <- regsubsets(x=form,data=df)
summary(mregsub)
summary(mregsub)$which
summary(mregsub)$adjr2
summary(mregsub)$rsq
summary(mregsub)$cp
summary(mregsub)$obj
summary(mregsub)$bic

mregsub
regOpt <- which.min(summary(mregsub)$bic)
regOpt <- 4
regOptActive <- names(coef(mregsub,regOpt)[-1])
form <- formula(paste(responsevar,"~",paste(regOptActive,collapse="+")))
# form <- formula("responsevar ~ P_Mex300em390 + P_Sag275and295 + P_Sag350and400 + P_HIOhno2008")
# form <- formula("responsevar ~ P_Mex300em390 + P_Sag275and295 + P_Sag350and400")
# form <- formula("responsevar ~ P_fDOMex370em460  + P_Sag275and295 + P_HIOhno2008")

# form <- formula("responsevar ~ P_fDOMex370em460  + P_Sag350and400 + P_HIOhno2008")
# form <- formula("responsevar ~ P_fDOMex370em460")

mOpt <- lm(formula=form,data=df)
lmActive <- names(summary(mOpt)$coef[,1])[-1]
adjR2Opt <- summary(mOpt)$adj.r.squared

predicted <- predict(mOpt)
##
##################################################

axis.limits<- range(c(observed,predicted))

plot(observed,predicted,ylim=axis.limits,xlim=axis.limits,col=df$plotColors,
     main="",pch=20,xlab="",ylab="")
abline(0,1,lty=2)
legend("topleft",legend=names(siteColors),col=siteColors,pch=20,cex=0.8)
mtext("Human Bacteroides Preliminary Regression with Optical Properties",side=3,line=3,font=2)
mtext(paste("adjR^2 =",round(adjR2Opt,2)),side=3,line=2,font=1)
mtext("Log10(Human Bacteroides predicted)",side=2,line=3,font=2)
mtext("Log10(Human Bacteroides observed)",side=1,line=3,font=2)
parmnames <- sub("P_","",lmActive)
mtext(paste(parmnames,collapse="; "),line=0.5)

################################################################################