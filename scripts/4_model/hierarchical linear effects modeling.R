# Reference to Sam's dissertation:  
# https://github.com/limnoliver/CSI-Nutrient-Time-Series/blob/master/Code/05_analysis_hlm.R

#Intro to Linear mixed models: https://stats.idre.ucla.edu/other/mult-pkg/introduction-to-linear-mixed-models/

library(lme4)
library(smwrBase)

#set data directories
raw.path <- "raw_data"
cached.path <- "cached_data"
summary.path <- "SummaryVariables"
summary.save <- "1_SummaryVariables"
cached.save <- "0_munge"

df <- readRDS(file.path(cached.path,summary.save,"dfOptP3P4CombinedContWQ.rds"))

response <- "logLachno"
which(substr(names(df),1,1)=="A")
AbsVars <- names(df)[c(61:138,232:240)]
FlVars <- names(df)[c(17:60,139:231)]
IVs <- c(AbsVars,FlVars)
df$sinDate <- fourier(df$psdate)[,1]
df$cosDate <- fourier(df$psdate)[,2]
df$UW <- ifelse(df$abbrev=="UW",1,0)
df$MC <- ifelse(df$abbrev=="MC",1,0)


df$season <- as.factor(floor((as.POSIXlt(df$psdate)$mon+1)/3))
df$season <- ifelse(as.POSIXlt(df$psdate)$mon == 11,1,df$season)

test <- df[,c("season","psdate")]
df$season <- ifelse(as.POSIXlt(df$psdate)$mon == 11,1,df$season)

df$CSOVol <- 0
df[which(df$GRnumber=="gr15841"),"CSOVol"] <- 594.7
df[which(df$GRnumber=="gr15887"),"CSOVol"] <- 524.9
df[which(df$GRnumber=="gr17417"),"CSOVol"] <- 341.2
df$CSO <- ifelse(df$CSOVol > 0,1,0)

IVs <- c("rA350_A370","rMrange.25_H2","B", "UW", "MC")

selectedRows <- which(df$abbrev %in% c("HW","MC","MW","UW"))
selectedRows <- which(df$abbrev %in% sites)
selectedRows <- which(df$abbrev %in% c("HW","MC","MW","UW","MF","CG","LD"))
selectedRows <- which(df$abbrev %in% c("MC","MW","UW"))


form <- paste(response," ~ ", 
paste(IVs,collapse = " + "),"+"," (",paste(IVs[1:3], collapse = " + "),"|","season",")")

# lmer(form,data=df)
# 
# lmer(lachno2 ~ rA350_A370 + B (rA350_A370 + B | season), data = df)

df$logLachno <- log10(df$lachno)
m <- lmer(logLachno ~ rA350_A370 + B + rMrange.25_H2 + UW + MC + (rA350_A370 + B + rMrange.25_H2 | season),data=df[selectedRows,])

m <- lmer(logLachno ~ T*sinDate + T*cosDate + F*cosDate + F*sinDate + Turbidity_mean + CSO +( T + F + Turbidity_mean | abbrev),data=df[selectedRows,])

# m <- lmer(logLachno ~ T + F + Turbidity_mean + CSO + ( T + F + Turbidity_mean | abbrev) + ( T + F + Turbidity_mean | season),data=df[selectedRows,])

predict(m)
summary(m)

#plot(fitted(m) ~ df[selectedRows,"logLachno"])

#Plot Observed vs Predicted from cv-determined model

colorOptions <- c("orange","yellow2","skyblue","black","springgreen4","blue","grey","darkorchid1")
names(colorOptions) <- unique(df$abbrev)

plotColors <- colorOptions[df[selectedRows,"abbrev"]]
levels(as.factor(df[selectedRows,"abbrev"]))

axis.log <- ""
axis.limits <- c(1,7)

par(mfcol=c(1,1))
plot(df[selectedRows,response],fitted(m),
     xlab="Observed",ylab="Predicted",col=plotColors,pch=20,log=axis.log,ylim=axis.limits,xlim=axis.limits)
abline(0,1)
#mtext(paste(Active.Coef.names[2:length(Active.Coef.names)],collapse=", "),side=3,line=1,cex=0.8)
mtext(paste("Linear mixed effects moedel for MMSD subwatersheds",response),side=3,line=2,font=2,cex=1)
legend(x="topleft",legend=names(colorOptions),col=colorOptions,pch=20,text.col=colorOptions,cex=0.7)

