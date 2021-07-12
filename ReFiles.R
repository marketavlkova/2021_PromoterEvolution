#!/usr/local/bin/Rscript
### function to order files acquired by Sys.glob function
### by numbers of wells instead of alphabetical order for MeanPlots.R script

reorder <- function(x, len) {
  ### check which column number has first data available
  name <- unlist(strsplit(x[1], split = ".", fixed = T))[1]
  well <- unlist(strsplit(name, "_"))[3]
  number <- unlist(strsplit(well, ""))[2:3]
  ### save well number
  if (number[1] == 0) {
    number <- number[2]
  } else {
    number <- as.numeric(number[1]) * 10 + as.numeric(number[2])
  }
  ### set value describing difference between column 1 and 1st column with data
  move <- as.numeric(number) - 1
  ### set output variable
  y <- character()
  for (i in 1:len) {
    for (xi in x) {
      ### extract important part from the filename
      name <- unlist(strsplit(xi, split = ".", fixed = T))[1]
      well <- unlist(strsplit(name, "_"))[3]
      number <- unlist(strsplit(well, ""))[2:3]
      ### save well number
      if (number[1] == 0) {
        number <- number[2]
      } else {
        number <- as.numeric(number[1]) * 10 + as.numeric(number[2])
      }
      ### when the well number correspond to "i"
      ### include it into output variable
      if ((as.numeric(number) - move) == i) {
        y <- c(y, xi)
      }
    }
  }
  return(y)
}
