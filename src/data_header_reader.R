#  data_header_reader.R
#  This file reads in .DAT and .CTL files created for ADMB .TPL files
#  The use must ensure not to have within line comments in the files (nope, fixed that)

# filename="PWS_ASA.dat"
# filename="PWS_ASA(age_comps).ctl"
# filename="PWS_ASA(disease).ctl"
# filename="PWS_ASA(surveys).ctl"
data_header_reader <- function(filename) {
  #  The user needs to make sure there is a blank line at the end of
  #  filename and that each data type (vector number or matrix) is
  #  separated by a blank line
  
  # This is kind of convoluted
  text <- readLines(filename)
  values <- grep("^[0-9]",text)
  signed.values <- grep("^[-]",text)
  read.these <- sort(c(values,signed.values))
  nlines <- length(text) 
  indices <- seq(1:nlines)
  indices <- indices[-read.these]
  first.differences <- c(diff(indices),5)# This accounts for the last data
  data.types <- length(first.differences[first.differences>1])
  
  headers <- list(text[indices],first.differences>1)
  return(headers)
}

