#!/usr/local/bin/Rscript
### function to extract fluorescence mean values from a set of FCS files
library("flowCore") ### handling FCS files
### incorporate gating function
if (!exists("foo", mode = "function")) source("FCgating.R")

### x is vector of FCS file names
### len is number of replicates per strain
exPeaks <- function(x, len) {
  ### set counter to track number of analyzed strains
  n <- 1
  y <- numeric()
  z <- numeric()
  ### loop through all files
  for (i in x) {
    ### read FCS file; tranformation argument is "linearize" by default
    ### to avoid conversion of some measured parameters into "a * 10^(x / R)"
    data <- read.FCS(i, alter.names = T)
    ### exclude wells with no cells
    all <- length(exprs(data$FSC.H))
    if (all <= 2500) {
      gfp.peak <- NA
    } else {
      ### gating by FCgating.R function
      cells <- FCgating(data)
      ### extract gfp values of cells in log10
      gfp.log <- log10(exprs(cells$GFP.H))
      ### replace infinite values by 0s
      gfp.log[is.infinite(gfp.log)] <- 0
      ### peak value of log10(gfp)
      denspl <- density(gfp.log)
      peak <- which.max(denspl$y)
      gfp.peak <- denspl$x[peak]
    }
    ### store the mean values in the output variable
    z <- c(z, gfp.peak)
    ### reset the counters when all values for the strain were plotted
    if (n < len) {
      n <- n + 1
    } else {
      n <- 1
    }
  }
  ### save all means as a matrix
  out <- matrix(z, ncol = len, byrow = T)
  colnames(out) <- c(1:len)
  rownames(out) <- c(1:(length(z)/len))
  write.csv(out, "AllPeaks.csv", row.names = T)
  ### count mean for each strain from all replicates
  n <- 1
  for (rep in z) {
    if (n == 1) {
      group <- rep
    } else {
      group <- c(group, rep)
    }
    if (n == len) {
      y <- c(y, mean(group, na.rm = T))
      n <- 1
    } else {
      n <- n + 1
    }
  }
  return(y)
}
