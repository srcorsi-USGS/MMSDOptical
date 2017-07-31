#Test regressions with common optical signals

library(glmnet)
library(smwrBase)
library(survival)
library(MASS)
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

#Define plotting colors
colorOptions <- c("orange","yellow2","skyblue","black","springgreen4","blue","grey","darkorchid1")
names(colorOptions) <- unique(df$abbrev)
selectedSiteColors <- colorOptions[selectedSites]
plotColors <- colorOptions[df[selectedRows,"abbrev"]]

#Define data set to use in modeling
df <- df[selectedRows,]

# Work with model #4 from HeuristicModels.R
# 
# 4. This model uses seasonality and water temp with peaks T and F plus CSO and site variables

#Original model:
m <- lm(log10(lachno2)~ Water_Temperature_mean  + CSO + T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW, data = df)
summary(m)
selectedRows <- c(1:dim(df)[1])
plotColors <- colorOptions[df[,"abbrev"]]

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)

#Examine errors
predicted <- predict(m,newdata = df)
RMSE <- mean(sqrt((df[,response] - predicted)^2));RMSE


#Split data into test and training sets and check RMSE
testRows <- sample(c(1:dim(df)[1]),size = round(dim(df)[1]*0.2))
dfTrain <- df[-testRows,]
dfTest <- df[testRows,]

m <- lm(log10(lachno2)~ Water_Temperature_mean  + CSO + T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW, data = dfTrain)
summary(m)

#Examine out-of-sample errors
predicted <- predict(m,newdata = dfTest)
RMSE <- mean(sqrt((df[testRows,response] - predicted)^2));RMSE

####
# Use caret to cross validate
set.seed(25)

#Fit model
form <- formula(loglachno2~ Water_Temperature_mean  + CSO + T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW)
m <- lm(form,data=df)
summary(m)
cM <- train(form, 
            data = df,
            method="lm",
            trControl = trainControl(
              method= "cv", number = 10,
              verboseIter = FALSE
            )
)

summary(cM)
print(cM)
cM

# Set up parameters
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)

#test with GBM also
cM <- train(loglachno2~ Water_Temperature_mean  + CSO + T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW, 
            data = df,
            method="gbm",
            trControl = fitControl,
            verbose = FALSE
)

summary(cM)
print(cM)
cM

#### 
# E coli models

#Fit model
cM <- train(logEcoli~ Water_Temperature_mean  + CSO + T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW, 
            data = df,
            method="lm",
            trControl = trainControl(
              method= "repeatedcv", number = 5,
              repeats = 10,verboseIter = FALSE
            )
)

summary(cM)
print(cM)
cM

###
# Examine variable importance

set.seed(7)
# load the library
library(mlbench)
library(caret)
# train the model
cM <- train(loglachno2~ Water_Temperature_mean  + CSO + T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW, 
            data = df,
            method="lm",
            trControl = trainControl(
              method= "repeatedcv", number = 10,
              repeats = 10,verboseIter = FALSE
            )
)

# estimate variable importance
importance <- varImp(cM, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(importance)


###
#Feature selection 
# ensure the results are repeatable
set.seed(7)

# define the control using a random forest selection function
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
form <- formula("loglachno2~ Water_Temperature_mean  + CSO + T*sinDate + T*cosDate + F*sinDate + F*cosDate + MC + UW")

x <- model.matrix(form,data=df)
results <- rfe(x[,-1], df[,response], sizes=c(3:13), rfeControl=control,metric="RMSE",maximize=FALSE)
# summarize the results
print(results)
# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))


######
### Cross validate final selected model

#Define variables to use in final model
selectedVars <- predictors(results)
#Add corresponding sin and cos terms
selectedVars <- c(selectedVars,"cosDate:F","T:sinDate", "MC")

form <- formula(paste("loglachno2 ~",paste(selectedVars,collapse = " + ")))
#Fit model
cM <- train(form, 
            data = df,
            method="lm",
            trControl = trainControl(
              method= "repeatedcv", number = 5,
              repeats = 10,verboseIter = FALSE
            )
)

summary(cM)
print(cM)
cM

m <- lm(form, data = df)

plotModel(m=m,df=df,response=response,selectedRows = c(1:dim(df)[1]),
          colorOptions = selectedSiteColors,plotColors = plotColors,
          pch=20)
summary(m)

