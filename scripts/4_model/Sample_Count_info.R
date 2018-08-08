
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

sum(!is.na(df$lachno2))
sum(!is.na(df$A254))
