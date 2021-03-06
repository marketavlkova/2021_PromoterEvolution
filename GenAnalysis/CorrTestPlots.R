#!/usr/local/bin/Rscript
### Script plotting Supp Figure 1a and Supp Figure 1b
### example of use: ./CorrTestPlots.R

library('scales')

### save output directory
dir <- dirname(getwd())

### load data file
data <- read.csv('output/DataIGRs&Genes.csv', header = T)
rownames(data) <- data[, 1]
data <- data[, -1]

### create list to customize promoter names for legend
prom.th <- list('aldAp' = 'aldA',
                'yhjXp' = 'yhjX',
                'mtrp2' = 'mtr',
                'aceBp' = 'aceB',
                'lacZp1' = 'lacZ',
                'dctAp' = 'dctA',
                'cddp' = 'cdd',
                'ptsGp2' = 'ptsG',
                'purAp' = 'purA',
                'tpiAp2' = 'tpiA')
### define colors for plotting
cols <- c('red', 'blue', 'green', 'darkviolet', 'darkorange', 'gold', 'saddlebrown', 'pink', 'grey', 'darkgreen')
names(cols) <- names(prom.th)

### plotting
pdf(file = paste0(dir, '/SupplementaryFigure_1a.pdf'), width = 5, height = 5)
par(las = 1, mar = c(5.1, 4.1, 2.1, 2.1))
  plot(data[, 'NormIgrTh'][which(data[, 1] >= 130)],
        data[, 'NormPromTh'][which(data[, 1] >= 130)],
        xlim = c(0, 0.08), ylim = c(0, 0.08),
        xlab = '', ylab = '', col = alpha('black', 0.25),
        bg = alpha('cyan', 0.25), pch = 21, cex = 0.9)
  for (pr in names(prom.th)) {
    points(data[pr, 'NormIgrTh'], data[pr, 'NormPromTh'], col = cols[pr], pch = 16)
  }
  abline(0, 1, col = 'red')
  c <- cor.test(data[, 'NormIgrTh'][which(data[, 1] >= 130)],
                data[, 'NormPromTh'][which(data[, 1] >= 130)],
                method = 'spearman', exact = F)
  mtext(text = paste0('rho = ', round(c$estimate, digits = 3)), side = 3, line = -1, cex = 0.9)
  mtext(text = paste0('p = ', signif(c$p.value, digits = 2)), side = 3, line = -2, cex = 0.9)
  title(xlab = expression(paste(theta, ' in intergenic regions')), line = 2.5)
  title(ylab = expression(paste(theta, ' in promoters')), line = 3)
  legend('topleft', legend = parse(text = sprintf('italic(%s)', prom.th)), col = cols,
          pch = 16, title = 'Promoter', cex = 0.85)
dev.off()

pdf(file = paste0(dir, '/SupplementaryFigure_1b.pdf'), width = 5, height = 5)
par(las = 1, mar = c(5.1, 4.1, 2.1, 2.1))
  plot(data[, 'IgrPi'][which(data[, 1] >= 130)],
        data[, 'PromPi'][which(data[, 1] >= 130)],
        xlim = c(0, 0.12), ylim = c(0, 0.12),
        xlab = '', ylab = '', col = alpha('black', 0.25),
        bg = alpha('cyan', 0.25), pch = 21, cex = 0.9)
  for (pr in names(prom.th)) {
    points(data[pr, 'IgrPi'], data[pr, 'PromPi'], col = cols[pr], pch = 16)
  }
  abline(0, 1, col = 'red')
  c <- cor.test(data[, 'IgrPi'][which(data[, 1] >= 130)],
                data[, 'PromPi'][which(data[, 1] >= 130)],
                method = 'spearman', exact = F)
  mtext(text = paste0('rho = ', round(c$estimate, digits = 3)), side = 3, line = -1, cex = 0.9)
  mtext(text = paste0('p = ', signif(c$p.value, digits = 2)), side = 3, line = -2, cex = 0.9)
  title(xlab = expression(paste(pi, ' in intergenic regions')), line = 2.5)
  title(ylab = expression(paste(pi, ' in promoters')), line = 3)
  legend('topleft', legend = parse(text = sprintf('italic(%s)', prom.th)), col = cols,
          pch = 16, title = 'Promoter', cex = 0.85)
dev.off()

c <- cor.test(data[, 'IgrPi'], data[, 'PromPi'], method = 'spearman', exact = F)
cat(paste('Correlation between Pi and Theta in IGRs:',
          round(c$estimate, digits = 4),
          'p-value:', signif(c$p.value, digits = 2)), '\n')
