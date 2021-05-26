#!/usr/local/bin/Rscript
### Script plotting Supp Figure 1c and Supp Figure 1d
### example of use: ./CorrTestPlots.R

library('scales')

### save output directory
dir <- dirname(getwd())

### load data files
ide <- read.csv('output/APIprom.csv', header = T)
seg <- read.csv('output/SegSitesProms.csv', header = T)
var <- read.csv('output/VariantProms.csv', header = T)
### create one data variable from loaded data files
proms <- matrix(, ncol = 5, nrow = length(seg[, 1]))
proms <- cbind(seg[, 3][order(seg[, 1])], seg[, 4][order(seg[, 1])], var[, 2], (100 - ide[, 2]), (seg[, 2][order(seg[, 1])] / seg[, 4][order(seg[, 1])]))
rownames(proms) <- var[, 1]
colnames(proms) <- c('NoStrains', 'AlignLength', 'NoVariants', '100-API', 'PSS')

### load data files
ide <- read.csv('output/APIigr.csv', header = T)
seg <- read.csv('output/SegSitesIGRs.csv', header = T)
var <- read.csv('output/VariantIGRs.csv', header = T)
### create one data variable from loaded data files
igrs <- matrix(, ncol = 5, nrow = length(seg[, 1]))
igrs <- cbind(seg[, 3][order(seg[, 1])], seg[, 4][order(seg[, 1])], var[, 2], (100 - ide[, 2]), (seg[, 2][order(seg[, 1])] / seg[, 4][order(seg[, 1])]))
rownames(igrs) <- var[, 1]
colnames(igrs) <- c('NoStrains', 'AlignLength', 'NoVariants', '100-API', 'PSS')

### load data files
ide <- read.csv('output/APIgene.csv', header = T)
seg <- read.csv('output/SegSitesGenes.csv', header = T)
var <- read.csv('output/VariantGenes.csv', header = T)
### create one data variable from loaded data files
genes <- matrix(, ncol = 5, nrow = length(seg[, 1]))
genes <- cbind(seg[, 3][order(seg[, 1])], seg[, 4][order(seg[, 1])], var[, 2], (100 - ide[, 2]), (seg[, 2][order(seg[, 1])] / seg[, 4][order(seg[, 1])]))
rownames(genes) <- var[, 1]
colnames(genes) <- c('NoStrains', 'AlignLength', 'NoVariants', '100-API', 'PSS')

### find promoters have corresponding gene pulled out
both_names <- c()
data <- matrix(, nrow = 0, ncol = 7)
for (pr in rownames(proms)) {
  for (ge in rownames(genes)) {
    if (pr == ge) {
      both_names <- c(both_names, pr)
      data <- rbind(data, c(proms[pr, 1], proms[pr, 3], proms[pr, 4], proms[pr, 5], genes[pr, 1], genes[pr, 4], genes[pr, 5]))
    }
  }
}
rownames(data) <- both_names
colnames(data) <- c('PromStr', 'PromVers', 'PromAPI', 'PromPSS', 'GeneStr', 'GeneAPI', 'GenePSS')
### save this as csv file
write.csv(data, file = 'output/DataProm&Genes.csv', row.names = T)

### create table with promoters that have both promoter and gene in at least 130 strains
prs <- c()
common <- matrix(, nrow = 0, ncol = 7)
for (pr in 1:length(data[, 1])) {
  if (data[pr, 1] >= 130 && data[pr, 5] >= 130) {
    prs <- c(prs, rownames(data)[pr])
    common <- rbind(common, data[pr,])
  }
}
rownames(common) <- prs
colnames(common) <- colnames(data)
### save this as csv file
write.csv(common, file = 'output/DataProm&GenesCommon.csv', row.names = T)

### create list to customize promoter names for legend
prom.pss <- list('aldAp' = 'aldA',
                'yhjXp' = 'yhjX',
                'lacZp1' = 'lacZ',
                'aceBp' = 'aceB',
                'mtrp2' = 'mtr',
                'cddp' = 'cdd',
                'dctAp' = 'dctA',
                'ptsGp2' = 'ptsG',
                'purAp' = 'purA',
                'tpiAp2' = 'tpiA')
### define colors for plotting
cols <- c('red', 'blue', 'green', 'darkviolet', 'darkorange', 'gold', 'saddlebrown', 'pink', 'grey', 'darkgreen')
names(cols) <- names(prom.pss)

### plotting
pdf(file = paste0(dir, '/SupplementaryFigure_1c.pdf'), width = 5, height = 5)
par(las = 1, mar = c(5.1, 4.1, 2.1, 2.1))
  plot(proms[, 5][which(proms[, 1] >= 130)],
        igrs[, 5][which(igrs[, 1] >= 130)],
        xlim = c(0, 0.415), ylim = c(0, 0.415),
        xlab = '', ylab = '', col = alpha('black', 0.25), pch = 16)
  for (pr in names(prom.pss)) {
    points(proms[pr, 5], igrs[pr, 5], col = cols[pr], pch = 16)
  }
  abline(0, 1, col = 'red')
  title(xlab = 'PSS with 100bp ORFs', line = 2.5)
  title(ylab = 'PSS without 100bp ORFs', line = 2.5)
  legend('topleft', legend = parse(text = sprintf('italic(%s)', prom.pss)), col = cols,
          pch = 16, title = 'Promoter', cex = 0.85)
dev.off()

pdf(file = paste0(dir, '/SupplementaryFigure_1d.pdf'), width = 5, height = 5)
par(las = 1, mar = c(5.1, 4.1, 2.1, 2.1))
  plot(proms[, 4][which(proms[, 1] >= 130)],
        igrs[, 4][which(igrs[, 1] >= 130)],
        xlim = c(0, 12), ylim = c(0, 12),
        xlab = '', ylab = '', col = alpha('black', 0.25), pch = 16)
  for (pr in names(prom.pss)) {
    points(proms[pr, 4], igrs[pr, 4], col = cols[pr], pch = 16)
  }
  abline(0, 1, col = 'red')
  title(xlab = '100 - API with 100bp ORFs', line = 2.5)
  title(ylab = '100 - API without 100bp ORFs', line = 2.5)
  legend('topleft', legend = parse(text = sprintf('italic(%s)', prom.pss)), col = cols,
          pch = 16, title = 'Promoter', cex = 0.85)
dev.off()
