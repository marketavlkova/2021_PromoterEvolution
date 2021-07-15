#!/usr/local/bin/Rscript
### Script plotting Figure 1b and Supp Figure 1b
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
lib.pss <- list()
lib.api <- list()
for (bc in bc0) {
  pss <- c()
  api <- c()
  for (row in 1:length(clear[, 1])) {
    if (clear[row, 3] == bc) {
      pss <- c(pss, clear[row, 5])
      api <- c(api, clear[row, 4])
    }
  }
  lib.pss[[bc]] <- pss
  lib.api[[bc]] <- api
}

### compare each function group to the rest of the dataset
pvals.pss <- c()
pvals.api <- c()
for (l in names(lib.pss)) {
  query.pss <- lib.pss[[l]]
  query.api <- lib.api[[l]]
  rest.pss <- c()
  rest.api <- c()
  for (ll in names(lib.pss)) {
    if (l != ll) {
      rest.pss <- c(rest.pss, lib.pss[[ll]])
      rest.api <- c(rest.api, lib.api[[ll]])
    }
  }
  wt.pss <- wilcox.test(x = query.pss, y = rest.pss)
  pvals.pss <- c(pvals.pss, wt.pss$p.value)
  wt.api <- wilcox.test(x = query.api, y = rest.api)
  pvals.api <- c(pvals.api, wt.api$p.value)
}
names(pvals.pss) <- names(lib.pss)
names(pvals.api) <- names(lib.api)

### plotting
cols <- c('blue', 'orange', 'green', 'red', 'purple',
          'saddlebrown', 'pink', 'grey', 'yellowgreen',
          'cyan')
pdf(file = paste0(dir, '/Figure_1b.pdf'), width = 5, height = 5)
par(las = 1, mar = c(5.1, 7, 2.1, 2.1))
  for (l in 1:10) {
    len <- length(lib.pss[[l]])
    if (l == 1) {
      plot(lib.pss[[l]], jitter(rep(l, len), factor = (10 / l)),
          xlim = c(0, 0.45), ylim = c(0.75, 10.25),
          yaxt = 'n', xaxt = 'n', pch = 16, col = alpha(cols[l], 0.4),
          xlab = '', ylab = '', cex = 0.8)
    } else {
      points(lib.pss[[l]], jitter(rep(l, len), factor = (10 / l)),
            pch = 16, col = alpha(cols[l], 0.4), cex = 0.8)
    }
    if (pvals.pss[l] <= (0.05 / 10)) {
      text(x = 0.425, y = l, labels = signif(pvals.pss[l], digits = 2), font = 4)
    } else {
      text(x = 0.425, y = l, labels = signif(pvals.pss[l], digits = 2))
    }
    arrows(median(lib.pss[[l]]), l - 0.25, median(lib.pss[[l]]), l + 0.25, length = 0)
  }
  axis(side = 1, at = seq(0, 0.4, by = 0.1), labels = seq(0, 0.4, by = 0.1), cex.axis = 0.75)
  axis(side = 2, at = c(1:10), labels = names(lib.pss), cex.axis = 0.75)
  title(xlab = 'Proportion of segregating sites in IGRs', line = 2, cex.lab = 0.8)
dev.off()

pdf(file = paste0(dir, '/SupplementaryFigure_1b.pdf'), width = 5, height = 5)
par(las = 1, mar = c(5.1, 7, 2.1, 2.1))
  for (l in 1:10) {
    len <- length(lib.api[[l]])
    if (l == 1) {
      plot(lib.api[[l]], jitter(rep(l, len), factor = (10 / l)),
          xlim = c(0, 14), ylim = c(0.75, 10.25),
          yaxt = 'n', xaxt = 'n', pch = 16, col = alpha(cols[l], 0.4),
          xlab = '', ylab = '', cex = 0.8)
    } else {
      points(lib.api[[l]], jitter(rep(l, len), factor = (10 / l)),
            pch = 16, col = alpha(cols[l], 0.4), cex = 0.8)
    }
    if (pvals.api[l] <= (0.05 / 10)) {
      text(x = 13, y = l, labels = signif(pvals.api[l], digits = 2), font = 4)
    } else {
      text(x = 13, y = l, labels = signif(pvals.api[l], digits = 2))
    }
    arrows(median(lib.api[[l]]), l - 0.25, median(lib.api[[l]]), l + 0.25, length = 0)
  }
  axis(side = 1, at = seq(0, 12.5, by = 2), labels = seq(0, 12.5, by = 2), cex.axis = 0.75)
  axis(side = 2, at = c(1:10), labels = names(lib.api), cex.axis = 0.75)
  title(xlab = '100 - average pairwise identity in IGRs', line = 2, cex.lab = 0.8)
dev.off()
