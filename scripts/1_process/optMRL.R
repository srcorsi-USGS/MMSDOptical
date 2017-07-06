#' Define function to compute MRL
#' optMRL
#'
#' Determines minimum reporting levels (MRLs) for absorbance data using results 
#' of blanks samples: MRL = mean(blank results) + 3 * sd(blank results).
#'
#' @param df dataframe with wavelength as column 1, and 
#' absorbance coefficients in the remaining columns
#' @param Wavelength character string representing the name of the wavelength
#' column in dataframe df.
#' @param blankGRnums vector of lab ID numbers that represent the blank samples
#' to be used in determination of the MRLs
#' @examples 
#' @export 
#' @return 

optMRL <- function(df,Wavelength,blankGRnums) {
  # Generate data frame with information on the blank samples by wavelength and
  # compute the minimum reporting level based on mean + 3 * SD for the blank samples
  # If the mean is less than zero, set the MRL to 3* SD
  dfBlankSummary <- data.frame(Wavelength = df[,Wavelength])
  dfBlankSummary$mean <- apply(df[,blankGRnums],MARGIN = 1, mean, na.rm=TRUE)
  dfBlankSummary$max <- apply(df[,blankGRnums],MARGIN = 1, max, na.rm=TRUE)
  dfBlankSummary$min <- apply(df[,blankGRnums],MARGIN = 1, min, na.rm=TRUE)
  dfBlankSummary$sd <- apply(df[,blankGRnums],MARGIN = 1, sd, na.rm=TRUE)
  dfBlankSummary$MRL <- dfBlankSummary$mean + 3 * dfBlankSummary$sd
  dfBlankSummary$MRL <- ifelse(dfBlankSummary$mean < 0, 3 * dfBlankSummary$sd, dfBlankSummary$MRL)
  return(dfBlankSummary)
}


