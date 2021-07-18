#!/usr/local/bin/Rscript
### source: https://stackoverflow.com/questions/35194048/using-r-how-to-calculate-the-distance-from-one-point-to-a-line

dist2d <- function(x0, x1, x2) {
  v1 <- x1 - x2
  v2 <- x0 - x1
  m <- cbind(v1, v2)
  d <- abs(det(m)) / sqrt(sum(v1 * v1))
  return(d)
}
