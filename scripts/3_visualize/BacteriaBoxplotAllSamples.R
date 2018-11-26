
#Boxplots of human bacteria
###################


library(latticeExtra)

#source(paste(Rlocal,"/Scripts/fxn_drop.unused.factors.R",sep=""))

#set data directories
raw.path <- "raw_data"
cached.path <- "cached_data"
summary.path <- "SummaryVariables"
summary.save <- "1_SummaryVariables"
cached.save <- "0_munge"

df <- readRDS(file.path(cached.path,summary.save,"dfOptP3P4CombinedContWQ.rds"))
df.orig <- df

selectedSites <- c("MC","MW","UW")
selectedRows <- which(df$abbrev %in% selectedSites)
colorOptions <- c("orange","yellow2","skyblue","black","springgreen4","blue","grey","darkorchid1")
names(colorOptions) <- unique(df$abbrev)
selectedSiteColors <- colorOptions[selectedSites]

df <- df[selectedRows,]
plotColors <- colorOptions[df[selectedRows,"abbrev"]]

sites<- c("Underwood_Cr", "Menomonee_R_Wauwatosa", "Menomonee_R_Milwaukee")
names(sites) <- c("UW","MW","MC")
df$Site.Name <- sites[df$abbrev]
#unique(df$Site.abbreviation)
df$Site.abbreviation <- factor(df$abbrev,levels=c("UW", "MW", "MC"))

name.levels <- as.character()
for (i in 1:length(levels(df$Site.abbreviation))){
  name.current <- as.character(unique(df[which(df$Site.abbreviation==levels(df$Site.abbreviation)[i]),"Site"]))
  name.levels <- c(name.levels,name.current)
}
df$Site.Name <- factor(df$Site.Name,levels=sites)


unique(df[which(df$Site.abbreviation==levels(df$Site.abbreviation)[i]),"Site.Name"])

#for (i in 10:13)  df[,i] <-as.character(df[,i])

# df$BacHuman <- as.character(df$bacHum)
# df$BacHuman <- ifelse(df$BacHuman=="BLD",112,df$BacHuman)
# df$BacHuman <- ifelse(df$BacHuman==0,112,df$BacHuman)
# df$BacHuman <- ifelse(df$BacHuman==1,112,df$BacHuman)

df$BacHuman <- as.numeric(df$bacHum)

bwplot(df$BacHuman~df$Site.Name,
       scales=list(y=list(log=TRUE)))

levels(df$Site.abbreviation)


#ifelse(df$Sample.Type=="Event ","Event",df$Sample.Type)

# dfbf <- subset(df,Sample.Type=="Baseflow")
# dfevent <- subset(df,Sample.Type=="Event")
# 
# df2 <- rbind(dfbf,dfevent)
# df2 <- drop.unused.factors(factordf=df2)
# df2$Site.Name <- factor(df2$Site.Name,levels=names(sites))
# df2$Site.abbreviation <- factor(df2$Site.abbreviation,levels=c("MA", "PO", "CL","RO","RM","OC", "JI","EE"))
# 
# bwplot(df2$BacHuman~df2$Site.abbreviation|df2$Sample.Type,
#        scales=list(y=list(log=TRUE)))
# 
# ############ Generate boxplot by site #####################################################################


#df3 <- read.csv("M:/QW Monitoring Team/GLRI toxics/Data Analysis/Data/compiledData/bacteriaWideFull.csv")
siteOrder <- sites
df$Site.Name <- factor(df$Site.Name,levels=siteOrder)

par()$mar
par()$oma

pdf("cached_data/2_visualize/HumanBacBoxplot.pdf",width=11,height=8)
par(oma=c(8,0,0,0),mfrow=c(1,1))

bp <- boxplot(df$BacHuman~df$hydro_condition + df$Site.Name,col=c(colors()[32],8),log="y",las=2,xaxt="n")
mtext("Human-specific Bacteroides for 3 Menomonee River Sites",side=3,line=2,font=2)

axis(side=1,at=((1:3)*2-0.5),labels=levels(df$Site.Name),las=2)
legend(x="topleft",legend=c("Baseflow","Event"),col=c(colors()[32],8),pch=c(15,15),bty="n",inset=0.01)
mtext(paste( bp$n, sep = ""), at = seq_along(bp$n), line = 0.1, side = 3)
mtext("n = ",side=3,line=0.1,adj=0)
dev.off()
shell.exec("HumanBacBoxplot.pdf")
############################################################################################################
# E. coli culture


par(oma=c(5,0,0,0),mfrow=c(1,1))

bp <- boxplot((df3$valueToUse_EColi+1)~df3$flowCondition+df3$labelName,col=c(colors()[32],8),
              log="y",las=2,xaxt="n")
mtext("E. coli (culture) for 8 Great Lakes Tributaries",side=3,line=2,font=2)
abline(v=c(4.5,10.5),lty=2)
mtext(c("Ohio","Michigan","Wisconsin"),side=1,line=6,at=c(2,7.5,13.5),font=2)

axis(side=1,at=((1:8)*2-0.5),labels=levels(df3$labelName),las=2)
legend(x=10.5,y=13000,legend=levels(df3$Sample.Type),col=c(colors()[32],8),pch=c(15,15),bty="n",inset=0.01)
mtext(paste( bp$n, sep = ""), at = seq_along(bp$n), line = 0.1, side = 3)
mtext("n = ",side=3,line=0.1,adj=0)


############################################################################################################

df4 <- read.csv("M:/QW Monitoring Team/GLRI toxics/Data Analysis/Data/compiledData/mergeVirusBac.csv")
siteOrder <- c("Maumee","Portage","Clinton","Rouge","Raisin","Menominee","Manitowoc","Milwaukee")
df4$labelName <- factor(df4$labelName,levels=siteOrder)
begin <- which(names(df4)=="Adenovirus.C.D.F")
end <- which(names(df4)=="HumanRotavirus")
df4$HumanSum <- apply(df4[,begin:end],1,sum)

begin <- which(names(df4)=="BVDV1" )
end <- which(names(df4)=="MycobacteriumAvium")
df4$BovineSum <- apply(df4[,begin:end],1,sum)

par(oma=c(3,3,3,1),mar=c(1,4,2,0),mfrow=c(2,2))
labelsize<-0.9
df4$HumanOccur <- factor(ifelse(df4$HumanSum>0,"Detect","Nondetect"),levels=c("Detect","Nondetect"))
bp <- boxplot((df4$P_BacteroidesHuman+1)~df4$HumanOccur,col=c(colors()[32],8),
              log="y",las=2,xaxt="n")
mtext("Human Bacteroides",side=3,line=0.5,col="blue",font=2,cex=labelsize)
mtext(paste( bp$n, sep = ""), at = seq_along(bp$n), line = 0.1, side = 1)
mtext("n = ",side=1,line=0.1,adj=0)



df4$HumanOccur <- factor(ifelse(df4$HumanSum>0,"Detect","Nondetect"),levels=c("Detect","Nondetect"))
bp <- boxplot((df4$P_EnterococcusQPCR+1)~df4$HumanOccur,col=c(colors()[32],8),
              log="y",las=2,xaxt="n")
mtext("Enterococcus QPCR",side=3,line=0.5,col="blue",font=2,cex=labelsize)
legend("bottom",legend=levels(df4$HumanOccur),col=c(colors()[32],8),pch=c(15,15),bty="n",inset=0.01)
mtext(paste( bp$n, sep = ""), at = seq_along(bp$n), line = 0.1, side = 1)
mtext("n = ",side=1,line=0.1,adj=0)



df4$HumanOccur <- factor(ifelse(df4$HumanSum>0,"Detect","Nondetect"),levels=c("Detect","Nondetect"))
bp <- boxplot((df4$P_EColi+1)~df4$HumanOccur,col=c(colors()[32],8),
              log="y",las=2,xaxt="n")
mtext("E. coli culture",side=3,line=0.5,col="blue",font=2,cex=labelsize)
mtext(paste( bp$n, sep = ""), at = seq_along(bp$n), line = 0.1, side = 1)
mtext("n = ",side=1,line=0.1,adj=0)
axis(side=1,at=1:2,labels=levels(df4$HumanOccur),line=1,col=0)


df4$HumanOccur <- factor(ifelse(df4$HumanSum>0,"Detect","Nondetect"),levels=c("Detect","Nondetect"))
bp <- boxplot((df4$P_Enterococcus+1)~df4$HumanOccur,col=c(colors()[32],8),
              log="y",las=2,xaxt="n")
mtext("Enterococcus culture",side=3,line=0.5,col="blue",font=2,cex=labelsize)
mtext(paste( bp$n, sep = ""), at = seq_along(bp$n), line = 0.1, side = 1)
mtext("n = ",side=1,line=0.1,adj=0)
axis(side=1,at=1:2,labels=levels(df4$HumanOccur),line=1,col=0)


mtext("Human virus detections compared to FIB",side=3,line=1.5,font=2,outer=T)

#save(df2,file="bachumanJune2012.RData")