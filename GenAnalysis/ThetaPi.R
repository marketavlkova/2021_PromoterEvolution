#!/usr/local/bin/Rscript
### Calculation of Pi & Watterson's theta
### example of use: ./ThetaPi.R

### save output directory
dir <- dirname(getwd())

### load promoter data files
ide <- read.csv('output/APIprom.csv', header = T)
seg <- read.csv('output/SegSitesProms.csv', header = T)
var <- read.csv('output/VariantProms.csv', header = T)
### create one data variable from loaded data files
proms <- matrix(, ncol = 5, nrow = length(seg[, 1]))
proms <- cbind(seg[, 3][order(seg[, 1])], seg[, 4][order(seg[, 1])], var[, 2], ide[, 2], seg[, 2][order(seg[, 1])])
rownames(proms) <- var[, 1]
colnames(proms) <- c('NoStrains', 'AlignLength', 'NoVariants', 'API', 'SegSites')

### load IGR data files
ide <- read.csv('output/APIigr.csv', header = T)
seg <- read.csv('output/SegSitesIGRs.csv', header = T)
var <- read.csv('output/VariantIGRs.csv', header = T)
### create one data variable from loaded data files
igrs <- matrix(, ncol = 5, nrow = length(seg[, 1]))
igrs <- cbind(seg[, 3][order(seg[, 1])], seg[, 4][order(seg[, 1])], var[, 2], ide[, 2], seg[, 2][order(seg[, 1])])
rownames(igrs) <- var[, 1]
colnames(igrs) <- c('NoStrains', 'AlignLength', 'NoVariants', 'API', 'SegSites')

### load gene data files
ide <- read.csv('output/APIgene.csv', header = T)
seg <- read.csv('output/SegSitesGenes.csv', header = T)
var <- read.csv('output/VariantGenes.csv', header = T)
### create one data variable from loaded data files
genes <- matrix(, ncol = 5, nrow = length(seg[, 1]))
genes <- cbind(seg[, 3][order(seg[, 1])], seg[, 4][order(seg[, 1])], var[, 2], ide[, 2], seg[, 2][order(seg[, 1])])
rownames(genes) <- var[, 1]
colnames(genes) <- c('NoStrains', 'AlignLength', 'NoVariants', 'API', 'SegSites')

### find promoters that have corresponding gene pulled out
both_names <- c()
data <- matrix(, nrow = 0, ncol = 13)
for (pr in rownames(igrs)) {
  for (ge in rownames(genes)) {
    if (pr == ge) {
      both_names <- c(both_names, pr)
      data <- rbind(data, c(igrs[pr, 1:5], genes[pr, 1:2], genes[pr, 4:5], proms[pr, 1:2], proms[pr, 4:5]))
    }
  }
}
rownames(data) <- both_names
colnames(data) <- c('IgrStr', 'IgrLen', 'IgrVers', 'IgrAPI', 'IgrSS', 'GeneStr', 'GeneLen', 'GeneAPI', 'GeneSS', 'PromStr', 'PromLen', 'PromAPI', 'PromSS')

### Theta calculation
ThetaI <- c()
ThetaP <- c()
ThetaG <- c()
### loop through all rows in data
for (row in 1:length(data[, 1])) {
  ### save number of strains & number of segregating sites
  ni <- data[row, 'IgrStr']
  np <- data[row, 'PromStr']
  ng <- data[row, 'GeneStr']
  Si <- data[row, 'IgrSS']
  Sp <- data[row, 'PromSS']
  Sg <- data[row, 'GeneSS']
  ### calculate harmonic number
  for (i in 1:(ni - 1)) {
    if (i == 1) {
      a1 <- 1/i
    } else {
      a1 <- a1 + (1/i)
    }
  }
  for (i in 1:(np - 1)) {
    if (i == 1) {
      a2 <- 1/i
    } else {
      a2 <- a2 + (1/i)
    }
  }
  for (i in 1:(ng - 1)) {
    if (i == 1) {
      a3 <- 1/i
    } else {
      a3 <- a3 + (1/i)
    }
  }
  ThetaI <- c(ThetaI, Si / a1 / data[row, 'IgrLen'])
  ThetaP <- c(ThetaP, Sp / a2 / data[row, 'PromLen'])
  ThetaG <- c(ThetaG, Sg / a3 / data[row, 'GeneLen'])
}

out <- cbind(data, ThetaI, ThetaG, ThetaP)
colnames(out) <- c(colnames(data), 'NormIgrTh', 'NormGeneTh', 'NormPromTh')

### Pi caulcation
igr.pi <- (100 - out[, 'IgrAPI']) / 100
prom.pi <- (100 - out[, 'PromAPI']) / 100
gene.pi <- (100 - out[, 'GeneAPI']) / 100
data <- cbind(out, as.vector(igr.pi), as.vector(gene.pi), as.vector(prom.pi))
colnames(data) <- c(colnames(out), 'IgrPi', 'GenePi', 'PromPi')

### save this as csv file
write.csv(data, file = 'output/DataIGRs&Genes.csv', row.names = T)

### create table with promoters that have both promoter and gene in at least 130 strains
prs <- c()
common <- matrix(, nrow = 0, ncol = length(data[1, ]))
for (pr in 1:length(data[, 1])) {
  if (data[pr, 'PromStr'] >= 130 && data[pr, 'GeneStr'] >= 130) {
    prs <- c(prs, rownames(data)[pr])
    common <- rbind(common, data[pr,])
  }
}
rownames(common) <- prs
colnames(common) <- colnames(data)

write.csv(common, file = 'output/DataIGRs&GenesCommon.csv')
