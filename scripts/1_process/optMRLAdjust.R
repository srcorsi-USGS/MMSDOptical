#' optMRLAdjust
#'
#' Define function to adjust absorbance data for minimum reporting levels (MRLs).
#'
#' @param df dataframe with wavelength as column 1, and 
#' absorbance coefficients in the remaining columns
#' @param dfMRLs dataframe with one column to define wavelength and one column 
#' to define MRLs
#' @param Wavelength character string representing the name of the wavelength
#' column in dataframes df and optMRL.
#' @param sampleGRnums vector of lab ID numbers that represent the environmental
#' samples
#' @param multiplier The value to multiply the MRL by to set abs values (typically
#' 1.0, 0.5, or 0). The default is 1.0 to set values below the MRL to the MRL.
#' @examples 
#' @export 
#' @return 

optMRLAdjust <- function(df,dfMRLs,Wavelength,sampleGRnums,multiplier=1.0) {  
  #Generate data frame with adjusted values based on the MRL. Generate a second 
  #dataframe with remark columns indicating values that are less than the MRL
  df2 <- data.frame(Wavelength = df[,Wavelength])
  dfRemarks <- data.frame(Wavelength = df[,Wavelength])
  for(colName in sampleGRnums){
    df2 <- cbind(df2,ifelse(df[,colName] < dfMRLs[,"MRL"],dfMRLs[,"MRL"]*multiplier,df[,colName]))
    dfRemarks <- cbind(dfRemarks,ifelse(df[,colName] < dfMRLs[,"MRL"],paste("<",dfMRLs[,"MRL"]),df[,colName]))
  }
  
  names(df2)[-1] <- sampleGRnums
  names(dfRemarks)[-1] <- sampleGRnums
  
  return(list(df2,dfRemarks))
}

