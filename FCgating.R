#!/usr/local/bin/Rscript
### gating function for FC data

FCgating <- function(data) {

  ### get max kernel density values from scatter values
  fsc.dens <- density(as.vector(exprs(data$FSC.H)))
  ssc.dens <- density(as.vector(exprs(data$SSC.H)))
  fypeak <- which.max(fsc.dens$y)
  sypeak <- which.max(ssc.dens$y)
  f.peak <- fsc.dens$x[fypeak]
  s.peak <- ssc.dens$x[sypeak]

  ### create rectangle gate at based on the kernel density
  pregate <- rectangleGate("FSC.H" = c(f.peak - 2000, f.peak + 15000),
                          "SSC.H" = c(s.peak - 2000, s.peak + 6000))
  filt <- filter(data, pregate)
  data <- Subset(data, filt)

  ### then create an elipse gate from cells filtered by rectangle gate
  FSC.H <- as.vector(exprs(data$FSC.H))
  SSC.H <- as.vector(exprs(data$SSC.H))
  SC.H <- cbind(FSC.H, SSC.H)
  cv <- cov(SC.H)                         ### covariance matrix needed for elipse gate definition
  mn <- c(median(FSC.H), median(SSC.H))   ### mean to define the elipse gate
  elgate <- ellipsoidGate(.gate = cv, mean = mn, distance = 1)
  filt <- filter(data, elgate)
  cells <- Subset(data, filt)

}
