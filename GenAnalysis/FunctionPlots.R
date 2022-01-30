#!/usr/local/bin/Rscript
### Script plotting Figure 1c and 1d
### example of use: ./FunctionPlots.R

library('scales')

### save output directory
dir <- dirname(getwd())

### create table for plotting using only the main function groups
if (!file.exists('output/MultiFunBC0-clean.csv')) {
  data <- read.csv('output/MultiFunPlotData.csv')

  data <- cbind(data[3:5], data[14:17])
  clear <- data[1,]
  for (row in 1:length(data[, 1])) {
    hit <- 0
    if (length(clear[, 1]) == 1) {
      a <- data[row, ]
      rownames(a) <- 1
      b <- clear[1, ]
      if (all(a == b)) {
        hit <- 1
      }
    } else {
      for (rowc in 1:length(clear[, 1])) {
        a <- data[row, ]
        rownames(a) <- 1
        b <- clear[rowc, ]
        rownames(b) <- 1
        if (all(a == b)) {
          hit <- 1
        }
      }
    }
    if (hit == 0) {
      clear <- rbind(clear, data[row, ])
    }
  }

  write.csv(clear, file = "output/MultiFunBC0-clean.csv")
} else {
  clear <- read.csv('output/MultiFunBC0-clean.csv', header = T)
  clear <- clear[, -1]
}

### get all possible function groups as a vector
bc0 <- c()
for (bc in clear[, 3]) {
  if (!is.element(bc, bc0)) {
    bc0 <- c(bc0, bc)
  }
}

### create a list item for each function group
lib.th <- list()
lib.pi <- list()
for (bc in bc0) {
  pss <- c()
  api <- c()
  for (row in 1:length(clear[, 1])) {
    if (clear[row, 3] == bc) {
      pss <- c(pss, clear[row, 5])
      api <- c(api, clear[row, 4])
    }
  }
  lib.th[[bc]] <- pss
  lib.pi[[bc]] <- api
}

### compare each function group to the rest of the dataset
pvals.th <- c()
pvals.pi <- c()
for (l in names(lib.th)) {
  query.th <- lib.th[[l]]
  query.pi <- lib.pi[[l]]
  rest.th <- c()
  rest.pi <- c()
  for (ll in names(lib.th)) {
    if (l != ll) {
      rest.th <- c(rest.th, lib.th[[ll]])
      rest.pi <- c(rest.pi, lib.pi[[ll]])
    }
  }
  wt.th <- wilcox.test(x = query.th, y = rest.th)
  pvals.th <- c(pvals.th, wt.th$p.value)
  wt.pi <- wilcox.test(x = query.pi, y = rest.pi)
  pvals.pi <- c(pvals.pi, wt.pi$p.value)
}
names(pvals.th) <- names(lib.th)
names(pvals.pi) <- names(lib.pi)

### plotting
cols <- c('blue', 'orange', 'green', 'red', 'purple',
          'saddlebrown', 'pink', 'grey', 'yellowgreen',
          'cyan')
pdf(file = paste0(dir, '/Figure_1c.pdf'), width = 5, height = 5)
par(las = 1, mar = c(5.1, 7, 2.1, 2.1))
  for (l in 1:10) {
    len <- length(lib.th[[l]])
    if (l == 1) {
      plot(lib.th[[l]], jitter(rep(l, len), factor = (10 / l)),
          xlim = c(0, 0.086), ylim = c(0.75, 10.25),
          yaxt = 'n', xaxt = 'n', pch = 16, col = alpha(cols[l], 0.4),
          xlab = '', ylab = '', cex = 0.8)
    } else {
      points(lib.th[[l]], jitter(rep(l, len), factor = (10 / l)),
            pch = 16, col = alpha(cols[l], 0.4), cex = 0.8)
    }
    if (pvals.th[l] <= (0.05 / 10)) {
      text(x = 0.08, y = l, labels = signif(pvals.th[l], digits = 2), font = 4)
    } else {
      text(x = 0.08, y = l, labels = signif(pvals.th[l], digits = 2))
    }
    arrows(median(lib.th[[l]]), l - 0.25, median(lib.th[[l]]), l + 0.25, length = 0)
  }
  axis(side = 1, at = seq(0, 0.086, by = 0.02), labels = seq(0, 0.086, by = 0.02), cex.axis = 0.75)
  axis(side = 2, at = c(1:10), labels = names(lib.th), cex.axis = 0.75)
  title(xlab = expression(paste(theta, ' in intergenic regions')), line = 2, cex.lab = 0.8)
dev.off()

pdf(file = paste0(dir, '/Figure_1d.pdf'), width = 5, height = 5)
par(las = 1, mar = c(5.1, 7, 2.1, 2.1))
  for (l in 1:10) {
    len <- length(lib.pi[[l]])
    if (l == 1) {
      plot(lib.pi[[l]], jitter(rep(l, len), factor = (10 / l)),
          xlim = c(0, 0.15), ylim = c(0.75, 10.25),
          yaxt = 'n', xaxt = 'n', pch = 16, col = alpha(cols[l], 0.4),
          xlab = '', ylab = '', cex = 0.8)
    } else {
      points(lib.pi[[l]], jitter(rep(l, len), factor = (10 / l)),
            pch = 16, col = alpha(cols[l], 0.4), cex = 0.8)
    }
    if (pvals.pi[l] <= (0.05 / 10)) {
      text(x = 0.14, y = l, labels = signif(pvals.pi[l], digits = 2), font = 4)
    } else {
      text(x = 0.14, y = l, labels = signif(pvals.pi[l], digits = 2))
    }
    arrows(median(lib.pi[[l]]), l - 0.25, median(lib.pi[[l]]), l + 0.25, length = 0)
  }
  axis(side = 1, at = seq(0, 0.15, by = 0.04), labels = seq(0, 0.15, by = 0.04), cex.axis = 0.75)
  axis(side = 2, at = c(1:10), labels = names(lib.pi), cex.axis = 0.75)
  title(xlab = expression(paste(pi, ' in intergenic regions')), line = 2, cex.lab = 0.8)
dev.off()
