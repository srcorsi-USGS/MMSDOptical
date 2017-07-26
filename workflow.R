#Workflow for processing and modeling optical data from MMSD Phase III & IV
#virus/bacteria/optical monitoring projects.
#Objective: develop regressions that can predict wastewater presence and 
#magnitude  from optical signals and make field-level sensors

#set data directories
raw.path <- "raw_data"
cached.path <- "cached_data"
summary.path <- "SummaryVariables"
summary.save <- "1_SummaryVariables"
cached.save <- "0_ProcessMDLs"
script.path <- "scripts"
munge.path <- "0_munge"
process.path <- "1_process"

#Munge
source(file.path(script.path,munge.path,"compileVirusBactOptTracking.R"))

#Process
##Determine censored values based on minimum reporting levels
source(file.path(script.path,process.path,"absMRLDetermination.R"))
source(file.path(script.path,process.path,"flMRLDetermination.R"))

##Determine summary optical variables
source(file.path(script.path,process.path,"getSummaryOpticalMMSDP3.R"))
source(file.path(script.path,process.path,"getSummaryOpticalMMSDP4.R"))

##Combine summary optical variables from P3 and P4 into one dataframe for modeling
source(file.path(script.path,process.path,"combineP3P4.R"))




