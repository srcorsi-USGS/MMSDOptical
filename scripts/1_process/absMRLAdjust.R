#' absMRLAdjust
#'
#' Define function to adjust absorbance data for minimum reporting levels (MRLs).
#'
#' @param dfabs dataframe with wavelength as column 1, and 
#' absorbance coefficients in the remaining columns
#' @param dfMRLs dataframe with one column to define wavelength and one column 
#' to define MRLs
#' @param Wavelength character string representing the name of the wavelength
#' column in dataframes dfabs and absMRL.
#' @param sampleGRnums vector of lab ID numbers that represent the environmental
#' samples
#' @param multiplier The value to multiply the MRL by to set abs values (typically
#' 1.0, 0.5, or 0). The default is 1.0 to set values below the MRL to the MRL.
#' @examples 
#' @export 
#' @return 

absMRLAdjust <- function(dfabs,dfMRLs,Wavelength,sampleGRnums,multiplier=1.0) {  
  #Generate data frame with adjusted values based on the MRL. Generate a second 
  #dataframe with remark columns indicating values that are less than the MRL
  dfabs2 <- data.frame(Wavelength = dfabs[,Wavelength])
  dfabsRemarks <- data.frame(Wavelength = dfabs[,Wavelength])
  for(colName in sampleGRnums){
    dfabs2 <- cbind(dfabs2,ifelse(dfabs[,colName] < dfMRLs[,"MRL"],dfMRLs[,"MRL"]*multiplier,dfabs[,colName]))
    dfabsRemarks <- cbind(dfabsRemarks,ifelse(dfabs[,colName] < dfMRLs[,"MRL"],paste("<",dfMRLs[,"MRL"]),dfabs[,colName]))
  }
  
  names(dfabs2)[-1] <- sampleGRnums
  names(dfabsRemarks)[-1] <- sampleGRnums
  
  return(list(dfabs2,dfabsRemarks))
}

