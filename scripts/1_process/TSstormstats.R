
TSstormstats <- function (df, date = "pdate", varname, dates, starttime = "Ebpdate", 
          endtime = "Eepdate", stats.return = c("mean"), subdfvar = "", 
          subdfvalue = "", subdatesvar = "", subdatesvalue = "", out.varname = "") 
{
  stats.names <- c("mean", "max", "min", "median", "sum", "sd", 
                   "maxdiff", "difference", "nearest", "nearprev")
  stats.get <- data.frame(row.names = 1)
  stats.get[, stats.names[1:length(stats.names)]] <- FALSE
  nstats <- length(stats.return)
  stats.get[, stats.return[1:nstats]] <- TRUE
  stats.return <- names(stats.get[which(stats.get[1, ] == TRUE)])
  if (subdfvar != "") 
    df <- subset(df, df[, subdfvar] == subdfvalue)
  if (subdatesvar != "") 
    dates <- subset(dates, dates[, subdatesvar] == subdatesvalue)
  maxrows <- nrow(dates)
  varstats = data.frame(row.names = 1:maxrows)
  resultname = vector(mode = "character")
  varcols <- which(names(df) == varname)
  for (k in 1:length(varcols)) {
    varname <- names(df)[varcols[k]]
    for (i in 1:maxrows) {
      subdata <- df[which(df[, date] >= dates[i, starttime] & 
                            df[, date] <= dates[i, endtime]), ]
      if(sum(!is.na(subdata[,varname]))>0){
      if (stats.get[, "mean"]) 
        varstats[i, "mean"] <- mean(subdata[, varname], 
                                    na.rm = T)
      if (stats.get[, "max"]) 
        varstats[i, "max"] <- max(subdata[, varname], 
                                  na.rm = T)
      if (stats.get[, "min"]) 
        varstats[i, "min"] <- min(subdata[, varname], 
                                  na.rm = T)
      if (stats.get[, "median"]) 
        varstats[i, "median"] <- median(subdata[, varname], 
                                        na.rm = T)
      if (stats.get[, "sum"]) 
        varstats[i, "sum"] <- sum(subdata[, varname], 
                                  na.rm = T)
      if (stats.get[, "sd"]) 
        varstats[i, "sd"] <- sd(subdata[, varname], na.rm = T)
      if (stats.get[, "maxdiff"]) 
        varstats[i, "maxdiff"] <- max(subdata[, varname], 
                                      na.rm = T) - min(subdata[, varname], na.rm = T)
      if (stats.get[, "difference"]) 
        varstats[i, "difference"] <- subdata[nrow(subdata), 
                                             varname] - subdata[1, varname]
      }
        
      statsnames <- stats.return
    }
    rm(resultname)
    resultname <- paste(out.varname, "_", stats.return, sep = "")
    names(varstats) <- resultname
    dates <- cbind(dates, varstats)
    varstats <- varstats[, -(1:length(varstats))]
  }
  return(dates)
}