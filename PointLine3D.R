#!/usr/local/bin/Rscript
### source: https://stackoverflow.com/questions/35194048/using-r-how-to-calculate-the-distance-from-one-point-to-a-line

dist3d <- function(x0, x1, x2) {
  v1 <- x1 - x2
  v2 <- x0 - x1
  v3 <- cross3d_prod(v1, v2)
  area <- sqrt(sum(v3 * v3)) / 2
  d <- 2 * area / sqrt(sum(v1 * v1))
  return(d)
}

cross3d_prod <- function(v1, v2) {
  v3 <- vector()
  v3[1] <- v1[2] * v2[3] - v1[3] * v2[2]
  v3[2] <- v1[3] * v2[1] - v1[1] * v2[3]
  v3[3] <- v1[1] * v2[2] - v1[2] * v2[1]
  return(v3)
}
