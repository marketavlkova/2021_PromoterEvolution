#!/usr/local/bin/Rscript
### Script producing plots from phenotypic analysis
### example to run: ./PlotPhenotype.R

library('flowCore') ### handling FCS files
library('scales') ### for function alpha
library('scatterplot3d') ### for 3D plot
library('tibble')
library('ggplot2')
library('ggridges')
### incorporate function extracting modal population expression values
if (!exists('foo', mode = 'function')) source('PeaksEx.R')
### incorporate function extracting coefficient of variation values
if (!exists('foo', mode = 'function')) source('CvEx.R')
### incorporate file reordering function
if (!exists("foo", mode="function")) source("ReFiles.R")
### incorporate function calculating distance of point from a line in 2D (plasticity)
if (!exists("foo", mode = "function")) source("PointLine2D.R")
### incorporate function calculating distance of point from a line in 3D (plasticity)
if (!exists("foo", mode = "function")) source("PointLine3D.R")

### save list of all possible conditions
### which is further used as a dictionary
cond.ls <- list('_glucose_' = 'Glucose',
                  '_pyruvic-acid_' = 'Pyruvic acid',
                  '_02-pyruvic-acid_' = 'Pyruvic acid (0.2%)',
                  '_L-malic-acid_' = 'L-malic acid',
                  '_fucose_' = 'Fucose',
                  '_glycerol_' = 'Glycerol',
                  '_glycerol+AMP_' = 'Glycerol + AMP',
                  '_glycerol+cytidine_' = 'Glycerol + Cyt',
                  '_galactose_' = 'Galactose',
                  '_lactose_' = 'Lactose',
                  '_glucose+tryptophan_' = 'Glucose + Try',
                  '_glucose+phenylalanine_' = 'Glucose + Phe',
                  '_D-mannose_' = 'D-mannose',
                  '_L-arabinose_' = 'L-arabionse')
### create a similar list as above with shortcuts
short.ls <- list('_glucose_' = 'Glu',
                  '_pyruvic-acid_' = 'Pyr',
                  '_02-pyruvic-acid_' = 'Pyr2',
                  '_L-malic-acid_' = 'Mal',
                  '_fucose_' = 'Fuc',
                  '_glycerol_' = 'Gly',
                  '_glycerol+AMP_' = 'Amp',
                  '_glycerol+cytidine_' = 'Cyt',
                  '_galactose_' = 'Gal',
                  '_lactose_' = 'Lac',
                  '_glucose+tryptophan_' = 'Try',
                  '_glucose+phenylalanine_' = 'Phe',
                  '_D-mannose_' = 'Man',
                  '_L-arabinose_' = 'Ara')
### create a list of promoters in order
### from highest Theta to the lowest Theta
prom.th <- list('AldA' = 'aldA',
                'YhjX' = 'yhjX',
                'Mtr' = 'mtr',
                'AceB' = 'aceB',
                'LacZ' = 'lacZ',
                'DctA' = 'dctA',
                'Cdd' = 'cdd',
                'PtsG' = 'ptsG',
                'PurA' = 'purA',
                'TpiA' = 'tpiA')

##################################################
######## DATA LOADING AND PREPROCESSING ##########
##################################################

### set working directories
root.path <- getwd()
### get promoter names and directories with FC datafiles
proms <- list.files(path = 'PhenAnalysis/')
prom.path <- vector()
data.path <- vector()
for (p in proms) {
  fc <- paste0(root.path, '/PhenAnalysis/', p, '/')
  prom.path <- c(prom.path, fc)
  dates <- list.files(path = paste0(fc, 'FC/'))
  for (d in dates) {
    set <- list.files(path = paste0(fc, 'FC/', d, '/'))
    data.path <- c(data.path, paste0(fc, 'FC/', d, '/', set, '/'))
  }
}
names(prom.path) <- proms

### save number of FCS datafiles for segregating variants
### (excluding lacZ which has a different layout)
nseg <- c(32, 36, 20, 16, 28, 15, 9, 8, 24)
nprom <- 1
nrep <- 1
### loop through all datasets and get max kernel densities
### from GFP channel & coefficients of variation
### (i.e. modal population expression & standard deviation / modal population expression)
ndir <- length(unlist(strsplit(data.path[1], split = '/', fixed = T)))
for (path in data.path) {
  setwd(path)
  if (!file.exists('Peaks.csv')) {
    pdir <- unlist(strsplit(path, split = '/', fixed = T))[ndir - 3]
    cat(paste('Analazying data from promoter', pdir, '\n'))
    ### extract info about last directory
    end.dir <- unlist(strsplit(path, split = '/', fixed = T))[ndir]
    ### change number of files in seg. promoter dir. accordingly
    if (!startsWith(end.dir, 'pxxx')) {
      nfiles <- 96
    } else {
      nfiles <- nseg[nprom]
      if (nrep < 3) {
        nrep <- nrep + 1
      } else {
        nrep <- 1
        nprom <- nprom + 1
      }
    }
    cond <- unlist(strsplit(end.dir, split = '_', fixed = T))[3]
    ### get all files names ending with .fsc
    Files <- Sys.glob('*.fcs')
    InFiles <- reorder(Files, nfiles)
    ### get modal population expression (max kernel density from GFP channel)
    Fpeaks <- exPeaks(InFiles, (length(Files) / nfiles))
    Fpeaks <- matrix(Fpeaks, nrow = nfiles, ncol = 1, byrow = T)
    rownames(Fpeaks) <- paste(1:nfiles)
    colnames(Fpeaks) <- cond
    write.csv(Fpeaks, 'Peaks.csv', row.names = T)
    ### get coefficient of variation, i.e., standard deviation
    ### from GFP channel / modal population expression
    Fcvs <- exCvs(InFiles, (length(Files) / nfiles))
    Fcvs <- matrix(Fcvs, nrow = nfiles, ncol = 1, byrow = T)
    rownames(Fcvs) <- paste(1:nfiles)
    colnames(Fcvs) <- cond
    write.csv(Fcvs, 'Cvs.csv', row.names = T)
  }
}

setwd(root.path)
### read peak value and SNPs info files for all promoters in glucose
mut.ls <- list()
seg.ls <- list()
offset.m <- vector()
offset.n <- vector()
offnames <- vector()
p <- 1
s <- 1
n <- 1
### change nseg variable to number of segregating variants instead
nseg <- c(18, 25, 12, 7, 19, 10, 4, 3, 16)
for (path in data.path) {
  ### extract info about last directory and growth condition
  end.dir <- unlist(strsplit(path, split = '/', fixed = T))[ndir]
  cond <- unlist(strsplit(end.dir, split = '_', fixed = T))[3]
  ### get info about type of the promoters from directory name
  prom <- unlist(strsplit(path, split = '/', fixed = T))[ndir - 3]
  ### if directory contains seg. promoter variants, save modal
  ### expression and coefficient of variation values in seg.ls list
  if (startsWith(end.dir, 'pxxx')) {
    data.m <- read.csv(paste0(path, 'Peaks.csv'), header = T)
    data.m <- data.m[, -1]
    data.m <- data.m[which(!is.na(data.m))]
    if (grepl('Mtr', path)) {
      data <- c(data.m[1:5], data.m[7:length(data.m)])
      data.m <- data
    }
    data.n <- read.csv(paste0(path, 'Cvs.csv'), header = T)
    data.n <- data.n[, -1]
    data.n <- data.n[which(!is.na(data.n))]
    if (grepl('Mtr', path)) {
      data <- c(data.n[1:5], data.n[7:length(data.n)])
      data.n <- data
    }
    df <- cbind(data.m[1:nseg[s]], data.n[1:nseg[s]])
    seg.ls[[paste0(prom, '_', cond, '_')]] <- df
    ### save values for offset calculation (in both
    ### expression and coefficient of variation)
    pos.seg.m <- data.m[nseg[s] + 1]
    neg.seg.m <- data.m[nseg[s] + 2]
    mg1655.seg.m <- data.m[nseg[s]]
    pos.seg.n <- data.n[nseg[s] + 1]
    neg.seg.n <- data.n[nseg[s] + 2]
    mg1655.seg.n <- data.n[nseg[s]]
    ### and save corresponding values from random
    ### mutagenesis in mut.ls list
    path.mut <- data.path[which((grepl(prom, data.path) &
                                  grepl(paste0('_', cond, '_'), data.path, fixed = T) &
                                  !grepl('pxxx', data.path)))]
    data.m <- read.csv(paste0(path.mut, 'Peaks.csv'), header = T)
    data.m <- data.m[, -1]
    data.n <- read.csv(paste0(path.mut, 'CVs.csv'), header = T)
    data.n <- data.n[, -1]
    df <- cbind(c(data.m[1:92], data.m[96]), c(data.n[1:92], data.n[96]))
    snp <- read.csv(paste0(prom.path[p], 'SimMatrixMutagenesis.csv'), header = T)
    snp <- snp[, 1:2]
    snp <- snp[-1,]
    snp <- rbind(snp, c('MG1655', 0))
    mut.ls[[paste0(prom, '_', cond, '_')]] <- cbind(snp, df)
    ### save values for offset
    pos.mut.m <- data.m[94]
    neg.mut.m <- data.m[95]
    mg1655.mut.m <- data.m[96]
    pos.mut.n <- data.n[94]
    neg.mut.n <- data.n[95]
    mg1655.mut.n <- data.n[96]
    ### calculate offset between segregating and random
    ### datasets (in both expression and coefficient of variation)
    off.m <- mean(c((pos.mut.m - pos.seg.m), (neg.mut.m - neg.seg.m), (mg1655.mut.m - mg1655.seg.m)), na.rm = T)
    offset.m <- c(offset.m, off.m)
    off.n <- mean(c((pos.mut.n - pos.seg.n), (neg.mut.n - neg.seg.n), (mg1655.mut.n - mg1655.seg.n)), na.rm = T)
    offset.n <- c(offset.n, off.n)
    offnames <- c(offnames, paste0(prom, '_', cond, '_'))
    if (n < 3) {
      n <- n + 1
    } else {
      s <- s + 1
      p <- p + 1
      n <- 1
    }
  ### if it is lacZ promoter directory do the same as above
  ### just the microplate layout is different from the others
  ### with zero offset as it is a single dataset
  } else if (startsWith(end.dir, 'placZ')) {
    data.m <- read.csv(paste0(path, 'Peaks.csv'), header = T)
    data.m <- data.m[, -1]
    data.n <- read.csv(paste0(path, 'Cvs.csv'), header = T)
    data.n <- data.n[, -1]
    mut.lacZ <- cbind(c(data.m[1:8], data.m[10:18], data.m[20:21], data.m[23:31], data.m[59]),
                c(data.n[1:8], data.n[10:18], data.n[20:21], data.n[23:31], data.n[59]))
    snp <- read.csv(paste0(prom.path[p], 'SimMatrixMutagenesis.csv'), header = T)
    snp <- snp[, 1:2]
    snp <- snp[-1,]
    names <- snp[, 1][which(snp[, 2] > 0 & snp[, 2] < 4)]
    names <- c(names, 'MG1655')
    vals <- snp[, 2][which(snp[, 2] > 0 & snp[, 2] < 4)]
    vals <- c(vals, 0)
    df.mut <- cbind(names, vals, mut.lacZ)
    mut.ls[[paste0(prom, '_', cond, '_')]] <- as.data.frame(df.mut)
    seg.ls[[paste0(prom, '_', cond, '_')]] <- cbind(c(data.m[41:48], data.m[33], data.m[49:59]),
                              c(data.n[41:48], data.n[33], data.n[49:59]))
    ### no offset needed as this is a single dataset
    offset.m <- c(offset.m, 0)
    offset.n <- c(offset.n, 0)
    offnames <- c(offnames, paste0(prom, '_', cond, '_'))
    if (n < 3) {
      n <- n + 1
    } else {
      p <- p + 1
      n <- 1
    }
  }
}
names(offset.m) <- offnames
names(offset.n) <- offnames

##################################################
################ PLOTTING BEGINS #################
##################################################

##################################################
############### CORRELATION PLOTS ################
############ FIGURE 2 AND ED FIGURE 3 ############
##################################################

### calculate madian value and standard deviation in modal population expression
### for all segregating variants of each promoter and condition
ranges <- list()
for (pr in proms) {
  args <- offnames[which(grepl(pr, offnames))]
  for (a in args) {
    seg.ok <- seg.ls[[a]][, 1]
    ### remove MG1655 variant of mtr promoter from calculation (SNP in GFP)
    if (grepl('Mtr', a)) {
      ok <- seg.ok[1:(length(seg.ok) - 1)]
      seg.ok <- ok
    }
    vals <- c(sd(seg.ok), median(seg.ok))
    names(vals) <- c('range', 'median')
    ranges[[a]] <- vals
  }
}

### set number of cloned (varC) and total number (varT)
### of segregating promoter variants
varC <- c(18, 25, 12, 7, 20, 20, 10, 4, 3, 16)
varT <- c(26, 32, 18, 7, 26, 25, 14, 6, 3, 22)
### save Theta & Pi values for IGRs
Th <- c(0.03128615749, 0.04649894747, 0.01402122108, 0.01602425267, 0.02988129083, 0.03669189673, 0.01549964916, 0.01238768076, 0.001703512842, 0.04397005734)
Pi <- c(0.0268, 0.0463, 0.0057, 0.0021, 0.0282, 0.0341, 0.0123, 0.01, 0.0005, 0.0555)
names(varC) <- proms
names(varT) <- proms
names(Th) <- proms
names(Pi) <- proms

### create tibble for ridge plots
ridge <- matrix(, ncol = 6, nrow = 402)
colnames(ridge) <- c('cond', 'varC', 'varT', 'Th', 'Pi', 'mode')
i <- 1
co <- c()
for (pr in proms) {
  args <- offnames[which(grepl(pr, offnames))]
  for (a in args) {
    seg.ok <- seg.ls[[a]][, 1]
    for (item in seg.ok) {
      b <- paste(pr, short.ls[[paste0('_', unlist(strsplit(a, split = '_'))[2], '_')]])
      ridge[i, ] <- c(b, varC[pr], varT[pr], Th[pr], Pi[pr], item)
      i <- i + 1
    }
    co <- c(co, b)
  }
}
cd <- as.character(ridge[, 1])
vc <- as.numeric(ridge[, 2])
vt <- as.numeric(ridge[, 3])
ps <- as.double(ridge[, 4])
ai <- as.double(ridge[, 5])
me <- as.double(ridge[, 6])

ord <- c()
for (i in 1:length(varT)) {
  ord <- c(ord, rep(varT[i], 3))
}

ridge <- tibble(cond = factor(cd, levels = co[order(ord)]),
                varC = vc, Variants = vt, Th = ps, Pi = ai, mode = me)

cat(paste('Producing Figure 2c\n'))
pdf(file = 'Figure_2c.pdf', width = 9, height = 7)
  par(las = 1, mar = c(5.1, 3.1, 2.1, 2.1))

  ggplot(ridge, aes(x = mode, y = cond, fill = Variants)) +
  geom_density_ridges(alpha = 0.9, stat = 'binline', bins = 60) +
  theme_ridges() +
  ylab('') +
  xlab('Modal expression (log10, a.u.)')

dev.off()

### create matrix combining all the values
### obtained above for each promoter and condition
mm <- matrix(, ncol = 6, nrow = 30)
mm2 <- matrix(, ncol = 6, nrow = 10)
i <- 1
i2 <- 1
for (pr in names(prom.th)) {
  row <- c(varC[pr], varT[pr], Th[pr], Pi[pr])
  args <- offnames[which(grepl(pr, offnames))]
  mm[i,] <- c(row, ranges[[args[1]]][1], ranges[[args[1]]][2])
  mm[i+1,] <- c(row, ranges[[args[2]]][1], ranges[[args[2]]][2])
  mm[i+2,] <- c(row, ranges[[args[3]]][1], ranges[[args[3]]][2])
  i <- i + 3
  args <- paste0(pr, '_glucose_')
  mm2[i2, ] <- c(row, ranges[[args]][1], ranges[[args]][2])
  i2 <- i2 + 1
}
colnames(mm) <- c('VariantsC', 'VariantsT', 'Th', 'Pi', 'SD', 'Mean')
colnames(mm2) <- colnames(mm)
rownames(mm) <- c(rep(proms[1], 3), rep(proms[2], 3),
        rep(proms[3], 3), rep(proms[4], 3),
        rep(proms[5], 3), rep(proms[6], 3),
        rep(proms[7], 3), rep(proms[8], 3),
        rep(proms[9], 3), rep(proms[10], 3))
rownames(mm2) <- names(prom.th)

### define colors for plotting and legend
cols <- c(rep('red', 3), rep('blue', 3), rep('green', 3),
        rep('darkviolet', 3), rep('darkorange', 3),
        rep('gold', 3), rep('saddlebrown', 3),
        rep('pink', 3), rep('grey', 3), rep('darkgreen', 3))
cols.leg <- c('red', 'blue', 'green', 'darkviolet', 'darkorange', 'gold', 'saddlebrown', 'pink', 'grey', 'darkgreen')
names(cols) <- names(prom.th)

### plot Figure 2a (correlation of Theta with stdev
### in modal expression of seg. variants)
cat(paste('Producing Figure 2a\n'))
pdf(file = 'Figure_2a.pdf', width = 5, height = 5)
  par(las = 1, mar = c(5.1, 4.1, 2.1, 2.1))

  plot(x = mm[, 3], y = mm[, 5], xlim = c(0, 0.05), ylim = c(0, 0.3),
        xlab = '', ylab = '', col = cols, pch = 16)
  points(x = mm2[, 3], y = mm2[, 5], pch = 21, col = 1)
  c <- cor.test(mm[, 3], mm[, 5], method = 'spearman', exact = F)
  c2 <- cor.test(mm2[, 3], mm2[, 5], method = 'spearman', exact = F)
  print(paste('glucose: rho =', round(c2$estimate, digits = 3), 'p =', signif(c2$p.value, digits = 2)))
  mtext(text = paste0('rho = ', round(c$estimate, digits = 3)), side = 3, line = -1, cex = 0.9)
  mtext(text = paste0('p = ', signif(c$p.value, digits = 2)), side = 3, line = -2, cex = 0.9)
  title(xlab = expression(paste(theta, ' in promoters')), line = 2.5)
  title(ylab = 'Stdev in modal expression', line = 3)
  legend('topleft', legend = parse(text = sprintf('italic(%s)', prom.th)), col = cols.leg,
          pch = 16, title = 'Promoter')

dev.off()

### plot Figure 2b (correlation of total seg. variant number
### with stdev in modal expression of seg. variants)
cat(paste('Producing Figure 2b\n'))
pdf(file = 'Figure_2b.pdf', width = 5, height = 5)
  par(las = 1, mar = c(5.1, 4.1, 2.1, 2.1))

  plot(x = mm[, 2], y = mm[, 5], xlim = c(3, 32), ylim = c(0, 0.3),
        xlab = '', ylab = '', col = cols, pch = 16)
  points(x = mm2[, 2], y = mm2[, 5], pch = 21, col = 1)
  c <- cor.test(mm[, 2], mm[, 5], method = 'spearman', exact = F)
  c2 <- cor.test(mm2[, 2], mm2[, 5], method = 'spearman', exact = F)
  print(paste('glucose: rho =', round(c2$estimate, digits = 3), 'p =', signif(c2$p.value, digits = 2)))
  mtext(text = paste0('rho = ', round(c$estimate, digits = 3)), side = 3, line = -1, cex = 0.9)
  mtext(text = paste0('p = ', signif(c$p.value, digits = 2)), side = 3, line = -2, cex = 0.9)
  title(xlab = 'Total number of segregating variants', line = 2.5)
  title(ylab = 'Stdev in modal expression', line = 3)
  legend('topleft', legend = parse(text = sprintf('italic(%s)', prom.th)), col = cols.leg,
          pch = 16, title = 'Promoter')

dev.off()

### plot Extended Data Figure 3a (correlation of Pi
### with stdev in modal expression of seg. variants)
cat(paste('Producing Extended Data Figure 3a\n'))
pdf(file = 'ED_Figure_3a.pdf', width = 5, height = 5)
  par(las = 1, mar = c(5.1, 4.1, 2.1, 2.1))

  plot(x = mm[, 4], y = mm[, 5], xlim = c(0, 0.06), ylim = c(0, 0.3),
        xlab = '', ylab = '', col = cols, pch = 16)
  c <- cor.test(mm[, 4], mm[, 5], method = 'spearman', exact = F)
  mtext(text = paste0('rho = ', round(c$estimate, digits = 3)), side = 3, line = -1, cex = 0.9)
  mtext(text = paste0('p = ', signif(c$p.value, digits = 2)), side = 3, line = -2, cex = 0.9)
  title(xlab = expression(paste(pi, ' in promoters')), line = 2.5)
  title(ylab = 'Stdev in modal expression', line = 3)
  legend('topleft', legend = parse(text = sprintf('italic(%s)', prom.th)), col = cols.leg,
          pch = 16, title = 'Promoter')

dev.off()

### plot Extended Data Figure 3b (correlation of number of cloned seg. variants
### with stdev in modal expression of seg. variants)
cat(paste('Producing Extended Data Figure 3b\n'))
pdf(file = 'ED_Figure_3b.pdf', width = 5, height = 5)
  par(las = 1, mar = c(5.1, 4.1, 2.1, 2.1))

  plot(x = mm[, 1], y = mm[, 5], xlim = c(3, 25), ylim = c(0, 0.3),
        xlab = '', ylab = '', col = cols, pch = 16)
  c <- cor.test(mm[, 1], mm[, 5], method = 'spearman', exact = F)
  mtext(text = paste0('rho = ', round(c$estimate, digits = 3)), side = 3, line = -1, cex = 0.9)
  mtext(text = paste0('p = ', signif(c$p.value, digits = 2)), side = 3, line = -2, cex = 0.9)
  title(xlab = 'Number of cloned segregating variants', line = 2.5)
  title(ylab = 'Stdev in modal expression', line = 3)
  legend('topleft', legend = parse(text = sprintf('italic(%s)', prom.th)), col = cols.leg,
          pch = 16, title = 'Promoter')

dev.off()

##################################################
################# MAPPING PLOTS ##################
#################### FIGURE 3 ####################
##################################################

### plot upper part of Figure 3a (mapping the effect of random SNPs
### relative to MG1655 variant to position within promoter sequence)
### and calculate differences caused by the SNPs in aceB promoter
### (save into the 'comp' variable)
cat(paste('Producing Figure 3a.1\n'))
comp <- list()
pdf(file = 'Figure_3a.1.pdf', width = 15, height = 7)
  par(mfcol = c(3, 5),
      las = 1)

  ### loop through first 5 promoters
  for (pr in names(prom.th)[1:5]) {
    ### extract info about all three environments
    args <- offnames[which(grepl(pr, offnames))]
    for (i in 1:3) {
      a <- args[i]
      arg <- paste0('_', unlist(strsplit(a, split = '_'))[2], '_')
      if (i == 1) {
        conds <- cond.ls[[which(endsWith(names(cond.ls), arg))]]
      } else {
        conds <- c(conds, cond.ls[[which(endsWith(names(cond.ls), arg))]])
      }
    }
    names(conds) <- args
    ### obtain data about SNP positions for each random variant
    snp_map <- read.csv(paste0(prom.path[pr], '1SNPmap.csv'), header = T)
    len <- length(snp_map[, 2])
    ### obtain info about promoter annotations (TF and so on)
    anns <- read.csv(paste0(prom.path[pr], 'AnnotationsBasic.csv'), header = T)
    ### loop through all three envrionments
    for (a in args) {
      ### save modal expression values for random variants and the MG1655 varaint
      pl <- as.numeric(mut.ls[[a]][, 3][which(as.numeric(mut.ls[[a]][, 2]) == 1)])
      mg <- as.numeric(mut.ls[[a]][length(mut.ls[[a]][, 3]), 3])
      ### in the case of aceB promoter calculate differences in
      ### modal expression relative to MG1655 variant for Figure 3b
      if (grepl('AceB', a)) {
        comp[[a]] <- (pl - mg)
      }
      ### plotting
      if (a == args[1]) {
        par(mar = c(0, 4, 3, 1))
        plot(3, 3, type = 'n', xlim = c(1, snp_map[len, 2]),
              ylim = c(2, 5), xaxt = 'n',
              main = parse(text = sprintf('italic(%s)', prom.th[pr])),
              cex.main = 1.5, xlab = '', ylab = '')
      } else if (a == args[2]) {
        par(mar = c(1.5, 4, 1.5, 1))
        plot(3, 3, type = 'n', xlim = c(1, snp_map[len, 2]),
              ylim = c(2, 5), xaxt = 'n',
              cex.main = 1.5, xlab = '', ylab = '')
      } else {
        par(mar = c(3, 4, 0, 1))
        plot(3, 3, type = 'n', xlim = c(1, snp_map[len, 2]),
              ylim = c(2, 5), xaxt = 'n',
              cex.main = 1.5, xlab = '', ylab = '')
      }
      for (an in 1:length(anns[, 1])) {
        rect(anns[an, 3], 2, anns[an, 4], 5, col = alpha(anns[an, 5], 0.25), border = NA)
      }
      mtext(conds[which(names(conds) == a)], side = 3, line = -1.5, adj = 0.01)
      arrows(1, mg, snp_map[len, 2], mg, length = 0)
      for (p in 1:length(pl)) {
        if (!is.na(pl[p])) {
          arrows(snp_map[p, 2], mg, snp_map[p, 2], pl[p],
                length = 0, col = alpha('black', 0.8), lwd = 1.5)
        }
      }
      start <- anns[length(anns[, 1]), 3]
      end <- anns[length(anns[, 1]), 4]
      positions <- seq(0 - floor(start / 50) * 50, floor((end - start) / 50) * 50, 50)
      if (a == args[3]) {
        axis(side = 1, at = seq(start - floor(start / 50) * 50, floor(end / 50) * 50 + 10 * floor((end - start) / 50) + 1, 50),
              labels = positions, cex.axis = 1)
      } else {
        axis(side = 1, at = seq(start - floor(start / 50) * 50, floor(end / 50) * 50 + 10 * floor((end - start) / 50) + 1, 50),
              labels = rep('', length(positions)), cex.axis = 1)
      }
      if (pr == names(prom.th)[1]) {
        title(ylab = 'Modal expression (log10, a.u.)', line = 2.5, cex.lab = 1.2)
      }
    }
  }

dev.off()

### plot upper part of Figure 3b (comparing differences in expression
### relative to MG1655 between conditions in aceB promoter)
cat(paste('Producing Figure 3b.1\n'))
pdf(file = 'Figure_3b.1.pdf', width = 9, height = 3)
  par(mfrow = c(1, 3),
      las = 1)

  all <- length(comp[[1]])
  plot(comp[[1]], comp[[2]], pch = 16,
        xlim = c(-2.25, 2.25), ylim = c(-2.25, 2.25),
        col = alpha('black', 0.3), xlab = '', ylab = '')
  title(xlab = 'Glucose', line = 2.5, cex.lab = 1.2)
  title(ylab = 'Pyruvic acid', line = 2.5, cex.lab = 1.2)
  abline(h = 0, col = 'blue', lty = 3)
  abline(v = 0, col = 'blue', lty = 3)
  both <- cbind(comp[[1]], comp[[2]])
  tl <- both[, 2][which(both[, 1] < 0)]
  tl <- tl[which(tl > 0)]
  tr <- both[, 2][which(both[, 1] > 0)]
  tr <- tr[which(tr > 0)]
  bl <- both[, 2][which(both[, 1] < 0)]
  bl <- bl[which(bl < 0)]
  br <- both[, 2][which(both[, 1] > 0)]
  br <- br[which(br < 0)]
  text(x = -1, y = 1, labels = round(length(tl) / all, digits = 2))
  text(x = 1, y = 1, labels = round(length(tr) / all, digits = 2))
  text(x = -1, y = -1, labels = round(length(bl) / all, digits = 2))
  text(x = 1, y = -1, labels = round(length(br) / all, digits = 2))

  plot(comp[[1]], comp[[3]], pch = 16,
        xlim = c(-2.25, 2.25), ylim = c(-2.25, 2.25),
        col = alpha('black', 0.3), xlab = '', ylab = '',
        main = parse(text = sprintf('italic(%s)', 'aceB')))
  title(xlab = 'Glucose', line = 2.5, cex.lab = 1.2)
  title(ylab = 'L-malic acid', line = 2.5, cex.lab = 1.2)
  abline(h = 0, col = 'blue', lty = 3)
  abline(v = 0, col = 'blue', lty = 3)
  both <- cbind(comp[[1]], comp[[3]])
  tl <- both[, 2][which(both[, 1] < 0)]
  tl <- tl[which(tl > 0)]
  tr <- both[, 2][which(both[, 1] > 0)]
  tr <- tr[which(tr > 0)]
  bl <- both[, 2][which(both[, 1] < 0)]
  bl <- bl[which(bl < 0)]
  br <- both[, 2][which(both[, 1] > 0)]
  br <- br[which(br < 0)]
  text(x = -1, y = 1, labels = round(length(tl) / all, digits = 2))
  text(x = 1, y = 1, labels = round(length(tr) / all, digits = 2))
  text(x = -1, y = -1, labels = round(length(bl) / all, digits = 2))
  text(x = 1, y = -1, labels = round(length(br) / all, digits = 2))

  plot(comp[[2]], comp[[3]], pch = 16,
        xlim = c(-2.25, 2.25), ylim = c(-2.25, 2.25),
        col = alpha('black', 0.3), xlab = '', ylab = '')
  title(xlab = 'Pyruvic acid', line = 2.5, cex.lab = 1.2)
  title(ylab = 'L-malic acid', line = 2.5, cex.lab = 1.2)
  abline(h = 0, col = 'blue', lty = 3)
  abline(v = 0, col = 'blue', lty = 3)
  both <- cbind(comp[[2]], comp[[3]])
  tl <- both[, 2][which(both[, 1] < 0)]
  tl <- tl[which(tl > 0)]
  tr <- both[, 2][which(both[, 1] > 0)]
  tr <- tr[which(tr > 0)]
  bl <- both[, 2][which(both[, 1] < 0)]
  bl <- bl[which(bl < 0)]
  br <- both[, 2][which(both[, 1] > 0)]
  br <- br[which(br < 0)]
  text(x = -1, y = 1, labels = round(length(tl) / all, digits = 2))
  text(x = 1, y = 1, labels = round(length(tr) / all, digits = 2))
  text(x = -1, y = -1, labels = round(length(bl) / all, digits = 2))
  text(x = 1, y = -1, labels = round(length(br) / all, digits = 2))

dev.off()

### plot lower part of Figure 3a (mapping the effect of random SNPs
### relative to MG1655 variant to position within promoter sequence)
### and calculate differences caused by the SNPs in dctA promoter
### (save into the 'comp' variable)
cat(paste('Producing Figure 3a.2\n'))
comp <- list()
pdf(file = 'Figure_3a.2.pdf', width = 15, height = 7)
  par(mfcol = c(3, 5),
      las = 1)

  ### loop through first 5 promoters
  for (pr in names(prom.th)[6:10]) {
    ### extract info about all three environments
    args <- offnames[which(grepl(pr, offnames))]
    for (i in 1:3) {
      a <- args[i]
      arg <- paste0('_', unlist(strsplit(a, split = '_'))[2], '_')
      if (i == 1) {
        conds <- cond.ls[[which(endsWith(names(cond.ls), arg))]]
      } else {
        conds <- c(conds, cond.ls[[which(endsWith(names(cond.ls), arg))]])
      }
    }
    names(conds) <- args
    ### obtain data about SNP positions for each random variant
    snp_map <- read.csv(paste0(prom.path[pr], '1SNPmap.csv'), header = T)
    len <- length(snp_map[, 2])
    ### obtain info about promoter annotations (TF and so on)
    anns <- read.csv(paste0(prom.path[pr], 'AnnotationsBasic.csv'), header = T)
    ### loop through all three envrionments
    for (a in args) {
      ### save modal expression values for random variants and the MG1655 varaint
      pl <- as.numeric(mut.ls[[a]][, 3][which(as.numeric(mut.ls[[a]][, 2]) == 1)])
      mg <- as.numeric(mut.ls[[a]][length(mut.ls[[a]][, 3]), 3])
      ### in the case of dctA promoter calculate differences in
      ### modal expression relative to MG1655 variant for Figure 3b
      if (grepl('DctA', a)) {
        comp[[a]] <- (pl - mg)
      }
      ### plotting
      if (a == args[1]) {
        par(mar = c(0, 4, 3, 1))
        plot(3, 3, type = 'n', xlim = c(1, snp_map[len, 2]),
              ylim = c(2, 5), xaxt = 'n',
              main = parse(text = sprintf('italic(%s)', prom.th[pr])),
              cex.main = 1.5, xlab = '', ylab = '')
      } else if (a == args[2]) {
        par(mar = c(1.5, 4, 1.5, 1))
        plot(3, 3, type = 'n', xlim = c(1, snp_map[len, 2]),
              ylim = c(2, 5), xaxt = 'n',
              cex.main = 1.5, xlab = '', ylab = '')
      } else {
        par(mar = c(3, 4, 0, 1))
        plot(3, 3, type = 'n', xlim = c(1, snp_map[len, 2]),
              ylim = c(2, 5), xaxt = 'n',
              cex.main = 1.5, xlab = '', ylab = '')
      }
      for (an in 1:length(anns[, 1])) {
        rect(anns[an, 3], 2, anns[an, 4], 5, col = alpha(anns[an, 5], 0.25), border = NA)
      }
      mtext(conds[which(names(conds) == a)], side = 3, line = -1.5, adj = 0.01)
      arrows(1, mg, snp_map[len, 2], mg, length = 0)
      for (p in 1:length(pl)) {
        if (!is.na(pl[p])) {
          arrows(snp_map[p, 2], mg, snp_map[p, 2], pl[p],
                length = 0, col = alpha('black', 0.8), lwd = 1.5)
        }
      }
      start <- anns[length(anns[, 1]), 3]
      end <- anns[length(anns[, 1]), 4]
      positions <- seq(0 - floor(start / 50) * 50, floor((end - start) / 50) * 50, 50)
      if (a == args[3]) {
        axis(side = 1, at = seq(start - floor(start / 50) * 50, floor(end / 50) * 50 + 10 * floor((end - start) / 50) + 1, 50),
              labels = positions, cex.axis = 1)
      } else {
        axis(side = 1, at = seq(start - floor(start / 50) * 50, floor(end / 50) * 50 + 10 * floor((end - start) / 50) + 1, 50),
              labels = rep('', length(positions)), cex.axis = 1)
      }
      if (pr == names(prom.th)[6]) {
        title(ylab = 'Modal expression (log10, a.u.)', line = 2.5, cex.lab = 1.2)
      }
    }
  }

dev.off()

### plot lower part of Figure 3b (comparing differences in expression
### relative to MG1655 between conditions in dctA promoter)
cat(paste('Producing Figure 3b.2\n'))
pdf(file = 'Figure_3b.2.pdf', width = 9, height = 3)
  par(mfrow = c(1, 3),
      las = 1)

  all <- length(comp[[1]])
  plot(comp[[1]], comp[[2]], pch = 16,
        xlim = c(-1, 1), ylim = c(-1, 1),
        col = alpha('black', 0.3), xlab = '', ylab = '')
  title(xlab = 'Glucose', line = 2.5, cex.lab = 1.2)
  title(ylab = 'Pyruvic acid', line = 2.5, cex.lab = 1.2)
  abline(h = 0, col = 'blue', lty = 3)
  abline(v = 0, col = 'blue', lty = 3)
  both <- cbind(comp[[1]], comp[[2]])
  tl <- both[, 2][which(both[, 1] < 0)]
  tl <- tl[which(tl > 0)]
  tr <- both[, 2][which(both[, 1] > 0)]
  tr <- tr[which(tr > 0)]
  bl <- both[, 2][which(both[, 1] < 0)]
  bl <- bl[which(bl < 0)]
  br <- both[, 2][which(both[, 1] > 0)]
  br <- br[which(br < 0)]
  text(x = -0.5, y = 0.5, labels = round(length(tl) / all, digits = 2))
  text(x = 0.5, y = 0.5, labels = round(length(tr) / all, digits = 2))
  text(x = -0.5, y = -0.5, labels = round(length(bl) / all, digits = 2))
  text(x = 0.5, y = -0.5, labels = round(length(br) / all, digits = 2))

  plot(comp[[1]], comp[[3]], pch = 16,
        xlim = c(-1, 1), ylim = c(-1, 1),
        col = alpha('black', 0.3), xlab = '', ylab = '',
        main = parse(text = sprintf('italic(%s)', 'dctA')))
  title(xlab = 'Glucose', line = 2.5, cex.lab = 1.2)
  title(ylab = 'L-malic acid', line = 2.5, cex.lab = 1.2)
  abline(h = 0, col = 'blue', lty = 3)
  abline(v = 0, col = 'blue', lty = 3)
  both <- cbind(comp[[1]], comp[[3]])
  tl <- both[, 2][which(both[, 1] < 0)]
  tl <- tl[which(tl > 0)]
  tr <- both[, 2][which(both[, 1] > 0)]
  tr <- tr[which(tr > 0)]
  bl <- both[, 2][which(both[, 1] < 0)]
  bl <- bl[which(bl < 0)]
  br <- both[, 2][which(both[, 1] > 0)]
  br <- br[which(br < 0)]
  text(x = -0.5, y = 0.5, labels = round(length(tl) / all, digits = 2))
  text(x = 0.5, y = 0.5, labels = round(length(tr) / all, digits = 2))
  text(x = -0.5, y = -0.5, labels = round(length(bl) / all, digits = 2))
  text(x = 0.5, y = -0.5, labels = round(length(br) / all, digits = 2))

  plot(comp[[2]], comp[[3]], pch = 16,
        xlim = c(-1, 1), ylim = c(-1, 1),
        col = alpha('black', 0.3), xlab = '', ylab = '')
  title(xlab = 'Pyruvic acid', line = 2.5, cex.lab = 1.2)
  title(ylab = 'L-malic acid', line = 2.5, cex.lab = 1.2)
  abline(h = 0, col = 'blue', lty = 3)
  abline(v = 0, col = 'blue', lty = 3)
  both <- cbind(comp[[2]], comp[[3]])
  tl <- both[, 2][which(both[, 1] < 0)]
  tl <- tl[which(tl > 0)]
  tr <- both[, 2][which(both[, 1] > 0)]
  tr <- tr[which(tr > 0)]
  bl <- both[, 2][which(both[, 1] < 0)]
  bl <- bl[which(bl < 0)]
  br <- both[, 2][which(both[, 1] > 0)]
  br <- br[which(br < 0)]
  text(x = -0.5, y = 0.5, labels = round(length(tl) / all, digits = 2))
  text(x = 0.5, y = 0.5, labels = round(length(tr) / all, digits = 2))
  text(x = -0.5, y = -0.5, labels = round(length(bl) / all, digits = 2))
  text(x = 0.5, y = -0.5, labels = round(length(br) / all, digits = 2))

dev.off()

### produce legend for Figure 3a (mapping the effect of random SNPs
### relative to MG1655 variant to position within promoter sequence)
cat(paste('Producing Figure 3a.3\n'))
pdf(file = 'Figure_3a.3.pdf', width = 15, height = 3)

  plot(NULL, axes = F, ann = F, xlim = c(0, 1), ylim = c(0, 1))
  legend('top', legend = c('open reading frame', '-35 or -10 element', 'inducer TF', 'repressor TF', 'dual TF', 'repressor srRNA'),
          pch = 15, col = c(alpha('yellow', 0.25), alpha('grey', 0.25), alpha('green', 0.25), alpha('red', 0.25), alpha('blue', 0.25), alpha('magenta', 0.25)),
          title = 'Annotations', bg = 'white', cex = 1.3, horiz = T)

dev.off()

### sorting effect sizes of changes in expression
### based on whether a SNP is inside (annTrue) or
### or outside (annFalse) TF or RNAP binding site
annTrue <- c()
annFalse <- c()
for (pr in names(prom.th)) {
  ### obtain data about SNP positions for each random variant
  snp_map <- read.csv(paste0(prom.path[pr], '1SNPmap.csv'), header = T)
  ### obtain info about promoter annotations (TF and so on)
  anns <- read.csv(paste0(prom.path[pr], 'AnnotationsBasic.csv'), header = T)
  ### extract info about all three environments
  args <- offnames[which(grepl(pr, offnames))]
  for (i in 1:3) {
    a <- args[i]
    arg <- paste0('_', unlist(strsplit(a, split = '_'))[2], '_')
    if (i == 1) {
      conds <- cond.ls[[which(endsWith(names(cond.ls), arg))]]
    } else {
      conds <- c(conds, cond.ls[[which(endsWith(names(cond.ls), arg))]])
    }
  }
  ### exclude Glucose + Phenylalanine environment due to
  ### a SNP in GFP gene in MG1655 variant of this promoter
  if (grepl('Mtr', pr)) {
    top <- 2
  } else {
    top <- 3
  }
  ### loop through all environments
  ### (except for the exclusion mentioned above)
  for (n in 1:top) {
    a <- args[n]
    ### save modal expression values for random variants and the MG1655 varaint
    pl <- as.numeric(mut.ls[[a]][, 3][which(as.numeric(mut.ls[[a]][, 2]) == 1)])
    mg <- as.numeric(mut.ls[[a]][length(mut.ls[[a]][, 3]), 3])
    ### calculate differences in log expression caused by individual SNPs
    sizes <- abs(pl - mg)
    ### loop through all random variants
    for (snp in 1:length(snp_map[, 2])) {
      hit <- 0
      for (ann in rownames(anns)) {
        if (grepl('TF', anns[ann, 2]) || grepl('sigma70', anns[ann, 2])) {
          if (snp_map[snp, 2] >= anns[ann, 3] && snp_map[snp, 2] <= anns[ann, 4]) {
            hit <- 1
          }
        }
      }
      ### if the variant has a SNP inside TF or RNAP binding site
      ###  assing the expression fold change value to 'annTrue' variable
      if (hit == 1) {
        annTrue <- c(annTrue, sizes[snp])
      ### if the SNP is outside of any TF or RNAP binding site
      ### assign the expression fold change value to 'annFalse' variable
      } else {
        annFalse <- c(annFalse, sizes[snp])
      }
    }
  }
}
### test for differences in the fold change sizes
test <- wilcox.test(x = annTrue, y = annFalse, alternative = 'greater')
cat(paste0('The fold-changes in expression due to single SNPs are larger inside RNAP and TF binding sites (',
          round(median(annTrue), digits = 2),
          ') than ouside of them (',
          round(median(annFalse, na.rm = T), digits = 2),
          ') with significance: p = ',
          signif(test$p.value, digits = 3), '\n'))

##################################################
################ VARIATION PLOTS #################
#################### FIGURE 4 ####################
##################################################

### produce Figure 4 (comparing span of expression values
### between segregating and random variants)
cat(paste('Producing Figure 4\n'))
pdf(file = 'Figure_4.pdf', width = 9, height = 7)
  par(mfrow = c(3, 4),
      las = 1)

  ### loop through all promoters
  for (pr in names(prom.th)) {
    if (pr == names(prom.th)[1] || pr == names(prom.th)[5] || pr == names(prom.th)[9]) {
      par(mar = c(3, 4, 2, 0))
    } else if (pr == names(prom.th)[2] || pr == names(prom.th)[6] || pr == names(prom.th)[10]) {
      par(mar = c(3, 3, 2, 1))
    } else if (pr == names(prom.th)[3] || pr == names(prom.th)[7]) {
      par(mar = c(3, 2, 2, 2))
    } else {
      par(mar = c(3, 1, 2, 3))
    }
    ### extract info about all three environments
    args <- offnames[which(grepl(pr, offnames))]
    for (i in 1:3) {
      a <- args[i]
      arg <- paste0('_', unlist(strsplit(a, split = '_'))[2], '_')
      if (i == 1) {
        conds <- short.ls[[which(endsWith(names(short.ls), arg))]]
      } else {
        conds <- c(conds, short.ls[[which(endsWith(names(short.ls), arg))]])
      }
    }
    ### loop through all environments
    for (n in 1:3) {
      a <- args[n]
      ### get modal expression values for random variants
      ### and use offset so they are comparable to seg. variants
      mall <- mut.ls[[a]][, 3]
      mall <- as.numeric(mall) - offset.m[a]
      ### get modal expression value for seg. variants
      sall <- seg.ls[[a]][, 1]
      ### remove MG1655 variant of mtr promoter from calculation (SNP in GFP)
      if (grepl('Mtr', a)) {
        ok <- sall[1:(length(sall) - 1)]
        sall <- ok
      }
      ### combine both variant groups to allow comparison of variances
      dat <- c(sall, mall)
      grp <- c(rep('A', length(sall)), rep('B', length(mall)))
      mat <- cbind(dat, grp)
      ### test for differences in variances (Fligner-Kileen test of homogeneity of variances)
      test <- fligner.test(x = dat, g = grp)
      ### plotting
      if (a == args[1]) {
        plot(jitter(rep(n + 0.2, length(mall)), factor = c(3:1)[n]), mall,
              main = parse(text = sprintf('italic(%s)', prom.th[pr])), col = alpha('red', 0.2),
              xlim = c(0.5, 3.5), ylim = c(2, 5), pch = 16,
              xlab = '', ylab = '', xaxt = 'n')
      } else {
        points(jitter(rep(n + 0.2, length(mall)), factor = c(3:1)[n]), mall,
              pch = 16, col = alpha('red', 0.2))
      }
      if (grepl('Mtr', a)) {
        points(jitter(rep(n - 0.2, length(sall)), factor = c(3:1)[n]), sall,
              pch = 16, col = alpha('black', 0.2))
      } else {
        points(jitter(rep(n - 0.2, (length(sall) - 1)), factor = c(3:1)[n]), sall[1:length(sall) - 1],
              pch = 16, col = alpha('black', 0.2))
        points(n - 0.2, sall[length(sall)],
              pch = 16, col = alpha('green', 0.75))
      }
      if (test$p.value <= (0.05 / 3)) {
        text(x = n, y = 4.9, labels = signif(test$p.value, digit = 2), font = 4)
      } else {
        text(x = n, y = 4.9, labels = signif(test$p.value, digit = 2))
      }
    }
    axis(side = 1, at = c(1:3), labels = conds)
    if (pr == names(prom.th)[1] || pr == names(prom.th)[5] || pr == names(prom.th)[9]) {
      title(ylab = 'Modal expression (log10, a.u.)', line = 2.5)
    }
  }
  plot(NULL, axes = F, ann = F, xlim = c(0, 1), ylim = c(0, 1))
  legend('top', legend = c('Mutagenesis', 'Segregating', 'MG1655'),
          pch = 16, col = c(alpha('red', 0.2), alpha('black', 0.2), alpha('green', 0.75)),
          title = 'Promoter set', cex = 1.5)

dev.off()

##################################################
################## PLASTICITY  ###################
############# FIGURE 5, ED FIGURE 4 ##############
##################################################

### save all modal expression values for aceB promoter
### (all environments, both seg. and random variants)
pr <- 'AceB'
args <- offnames[which(grepl(pr, offnames))]
mut1 <- mut.ls[[args[1]]][, 3]
mut1 <- as.numeric(mut1) - offset.m[args[1]]
mut2 <- mut.ls[[args[2]]][, 3]
mut2 <- as.numeric(mut2) - offset.m[args[2]]
mut3 <- mut.ls[[args[3]]][, 3]
mut3 <- as.numeric(mut3) - offset.m[args[3]]
seg1 <- seg.ls[[args[1]]][, 1]
seg2 <- seg.ls[[args[2]]][, 1]
seg3 <- seg.ls[[args[3]]][, 1]
nseg <- length(seg1)
for (i in 1:3) {
  a <- args[i]
  args[i] <- paste0('_', unlist(strsplit(a, split = '_'))[2], '_')
  if (i == 1) {
    conds <- cond.ls[[which(endsWith(names(cond.ls), args[i]))]]
  } else {
    conds <- c(conds, cond.ls[[which(endsWith(names(cond.ls), args[i]))]])
  }
}
### set values for ylim and xlim in plotting
lims <- c(2, 5)

### produce Figure 5a (comparing expression values of aceB
### promoter between pyruvic and L-malic acid)
cat(paste('Producing Figure 5a\n'))
pdf(file = "Figure_5a.pdf", width = 5, height = 5)
  par(las = 1)

  plot(mut2, mut3, pch = 16, col = alpha('red', 0.2),
      ylim = c(2, 5), xlim = c(2, 5), xlab = '', ylab = '')
  points(seg2[1:nseg - 1], seg3[1:nseg - 1],
        pch = 16, col = alpha('black', 0.2))
  points(seg2[nseg], seg3[nseg],
        pch = 16, col = alpha('green', 0.75))
  abline(0, 1, lty = 3, col = 'blue')
  mtext(parse(text = sprintf('italic(%s)', prom.th[pr])), side = 3, line = 0.5)
  title(xlab = paste(conds[2]), line = 2)
  title(ylab = paste(conds[3]), line = 2.5)

dev.off()

### produce Figure 5b (comparing expression values of aceB
### promoter among all three environments -> 3D)
cat(paste('Producing Figure 5b\n'))
pdf(file = "Figure_5b.pdf", width = 5, height = 5)
  par(las = 1)

  source(paste0(root.path, '/addgrids3d.r'))
  sp <- scatterplot3d(mut1, mut2, mut3,
        pch = '', grid = F, box = F,
        xlab = conds[1], xlim = lims,
        zlab = conds[2], ylim = lims,
        ylab = conds[3], zlim = lims,
        angle = 55)
  addgrids3d(lims, lims, lims, grid = c('xy', 'xz', 'yz'), angle = 55)
  sp$points(mut1, mut2, mut3,
        col = alpha("red", 0.2), pch = 16)
  sp$points3d(seg1[1:nseg - 1], seg2[1:nseg - 1], seg3[1:nseg - 1],
        col = alpha("black", 0.2), pch = 16)
  sp$points3d(seg1[nseg], seg2[nseg], seg3[nseg],
        col = alpha("green", 0.75), pch = 16)
  sp$points3d(x = lims, y = lims, z = lims, type = 'l', lty = 3, col = 'blue')
  mtext(parse(text = sprintf('italic(%s)', prom.th[pr])), side = 3, line = 0.5)

dev.off()

### produce Figure 5c (comparing plasticity in 3D
### of all promoters between seg. and random variants)
cat(paste('Producing Figure 5c\n'))
pdf(file = "Figure_5c.pdf", width = 9, height = 5)
  par(las = 1)

  ### define points in 3D on an isospline to use when
  ### calculating the distance of datapoints from the isospline
  x1 <- rep(2, 3)
  x2 <- rep(5, 3)

  ### loop through all promoters
  for (n in 1:length(names(prom.th))) {
    pr <- names(prom.th)[n]
    ### extract info about all three environments
    args <- offnames[which(grepl(pr, offnames))]
    a <- args[1]
    b <- args[2]
    c <- args[3]
    ### get values from random variants and
    ### correct them using offset values
    mut1 <- mut.ls[[a]][, 3]
    mut1 <- as.numeric(mut1) - offset.m[a]
    mut2 <- mut.ls[[b]][, 3]
    mut2 <- as.numeric(mut2) - offset.m[b]
    mut3 <- mut.ls[[c]][, 3]
    mut3 <- as.numeric(mut3) - offset.m[c]
    mut <- c()
    ### loop through each random variant and calculate
    ### its minimal distance from isospline (= plasticity)
    for (m in 1:length(mut1)) {
      x0 <- c(mut1[m], mut2[m], mut3[m])
      mut <- c(mut, dist3d(x0, x1, x2))
    }
    nmut <- length(mut)
    ### get values from seg. variants
    seg1 <- seg.ls[[a]][, 1]
    seg2 <- seg.ls[[b]][, 1]
    seg3 <- seg.ls[[c]][, 1]
    ### remove MG1655 variant of mtr promoter from calculation (SNP in GFP)
    if (grepl('Mtr', pr)) {
      ok <- seg1[1:(length(seg1) - 1)]
      seg1 <- ok
      ok <- seg2[1:(length(seg2) - 1)]
      seg2 <- ok
      ok <- seg3[1:(length(seg3) - 1)]
      seg3 <- ok
    }
    seg <- c()
    ### loop through each seg. variant and calculate
    ### its minimal distance from isospline (= plasticity)
    for(s in 1:length(seg1)) {
      x0 <- c(seg1[s], seg2[s], seg3[s])
      seg <- c(seg, dist3d(x0, x1, x2))
    }
    nseg <- length(seg)
    ### test for differences in median plasticity between seg. and random variants
    test <- wilcox.test(x = mut, y = seg, exact = F)
    ### plotting
    if (pr == names(prom.th)[1]) {
      plot(jitter(rep(n - 0.2, nseg - 1), factor = 4 / n), seg[1:nseg - 1],
              xaxt = 'n', ylim = c(0, 2), xlim = c(0.5, 10.5),
              col = alpha('black', 0.2), pch = 16, xlab = '', ylab = '')
    } else {
      points(jitter(rep(n - 0.2, nseg - 1), factor = 4 / n), seg[1:nseg - 1],
              col = alpha('black', 0.2), pch = 16)
    }
    points(jitter(rep(n + 0.2, nmut), factor = 4 / n), mut,
            col = alpha('red', 0.2), pch = 16)
    if (grepl('Mtr', pr)) {
      points(n - 0.2, seg[nseg],
              col = alpha('black', 0.2), pch = 16)
    } else {
      points(n - 0.2, seg[nseg],
              col = alpha('green', 0.75), pch = 16)
    }
    arrows(c(n - 0.3, n + 0.1), c(median(seg, na.rm = T), median(mut, na.rm = T)),
          c(n - 0.1, n + 0.3), c(median(seg, na.rm = T), median(mut, na.rm = T)), length = 0)
    if (test$p.value <= 0.05) {
      text(x = n, y = 1.95, labels = signif(test$p.value, digits = 2), font = 4)
    } else {
      text(x = n, y = 1.95, labels = signif(test$p.value, digits = 2))
    }
  }
  axis(side = 1, at = c(1:10), labels = parse(text = sprintf('italic(%s)', prom.th)))
  title(ylab = 'Plasticity', line = 2.5)
  legend('right', legend = c('Mutagenesis', 'Segregating', 'MG1655'),
          pch = 16, col = c(alpha('red', 0.2), alpha('black', 0.2), alpha('green', 0.75)),
          title = 'Promoter set')

dev.off()

### save all modal expression values for aceB promoter
### (all environments, both seg. and random variants)
pr <- 'DctA'
args <- offnames[which(grepl(pr, offnames))]
mut1 <- mut.ls[[args[1]]][, 3]
mut1 <- as.numeric(mut1) - offset.m[args[1]]
mut2 <- mut.ls[[args[2]]][, 3]
mut2 <- as.numeric(mut2) - offset.m[args[2]]
mut3 <- mut.ls[[args[3]]][, 3]
mut3 <- as.numeric(mut3) - offset.m[args[3]]
seg1 <- seg.ls[[args[1]]][, 1]
seg2 <- seg.ls[[args[2]]][, 1]
seg3 <- seg.ls[[args[3]]][, 1]
nseg <- length(seg1)
for (i in 1:3) {
  a <- args[i]
  args[i] <- paste0('_', unlist(strsplit(a, split = '_'))[2], '_')
  if (i == 1) {
    conds <- cond.ls[[which(endsWith(names(cond.ls), args[i]))]]
  } else {
    conds <- c(conds, cond.ls[[which(endsWith(names(cond.ls), args[i]))]])
  }
}

### produce Figure 5d (comparing expression values of dctA
### promoter between pyruvic and L-malic acid)
cat(paste('Producing Figure 5d\n'))
pdf(file = "Figure_5d.pdf", width = 5, height = 5)
  par(las = 1)

  plot(mut2, mut3, pch = 16, col = alpha('red', 0.2),
      ylim = lims, xlim = lims, xlab = '', ylab = '')
  points(seg2[1:nseg - 1], seg3[1:nseg - 1],
        pch = 16, col = alpha('black', 0.2))
  points(seg2[nseg], seg3[nseg],
        pch = 16, col = alpha('green', 0.75))
  abline(0, 1, lty = 3, col = 'blue')
  mtext(parse(text = sprintf('italic(%s)', prom.th[pr])), side = 3, line = 0.5)
  title(xlab = conds[2], line = 2)
  title(ylab = conds[3], line = 2.5)

dev.off()

### define points in 2D on an isospline to use when
### calculating the distance of datapoints from the isospline
x1 <- rep(2, 2)
x2 <- rep(5, 2)

### produce Figure 5e (comparing plasticity in 2D
### of dctA promoter between seg. and random variants)
cat(paste('Producing Figure 5e\n'))
pdf(file = "Figure_5e.pdf", width = 5, height = 5)
  par(las = 1)

  ### loop through each random variant and calculate
  ### its minimal distance from isospline (= plasticity)
  m1 <- c()
  for (m in 1:length(mut1)) {
    x0 <- c(mut1[m], mut2[m])
    m1 <- c(m1, dist2d(x0, x1, x2))
  }
  m2 <- c()
  for (m in 1:length(mut2)) {
    x0 <- c(mut2[m], mut3[m])
    m2 <- c(m2, dist2d(x0, x1, x2))
  }
  m3 <- c()
  for (m in 1:length(mut3)) {
    x0 <- c(mut3[m], mut1[m])
    m3 <- c(m3, dist2d(x0, x1, x2))
  }
  nmut <- length(m1)
  ### loop through each seg. variant and calculate
  ### its minimal distance from isospline (= plasticity)
  s1 <- c()
  for(s in 1:length(seg1)) {
    x0 <- c(seg1[s], seg2[s])
    s1 <- c(s1, dist2d(x0, x1, x2))
  }
  s2 <- c()
  for(s in 1:length(seg2)) {
    x0 <- c(seg2[s], seg3[s])
    s2 <- c(s2, dist2d(x0, x1, x2))
  }
  s3 <- c()
  for(s in 1:length(seg3)) {
    x0 <- c(seg3[s], seg1[s])
    s3 <- c(s3, dist2d(x0, x1, x2))
  }
  nseg <- length(s1)
  ### test for differences in 2D plasticity between seg. and random variants
  t1 <- wilcox.test(x = m1, y = s1, exact = F)
  t2 <- wilcox.test(x = m2, y = s2, exact = F)
  t3 <- wilcox.test(x = m3, y = s3, exact = F)
  t <- c(t1$p.value, t2$p.value, t3$p.value)
  ### plotting
  plot(jitter(rep(0.8, nseg - 1), factor = 3 / 1), s1[1:nseg - 1],
          xaxt = 'n', ylim = c(0, 2), xlim = c(0.5, 3.5),
          col = alpha('black', 0.2), pch = 16, xlab = '', ylab = '')
  points(jitter(rep(1.8, nseg - 1), factor = 3 / 2), s2[1:nseg - 1],
          col = alpha('black', 0.2), pch = 16)
  points(jitter(rep(2.8, nseg - 1), factor = 3 / 3), s3[1:nseg - 1],
          col = alpha('black', 0.2), pch = 16)
  points(jitter(rep(1.2, nmut), factor = 3 / 1), m1,
          col = alpha('red', 0.2), pch = 16)
  points(jitter(rep(2.2, nmut), factor = 3 / 2), m2,
          col = alpha('red', 0.2), pch = 16)
  points(jitter(rep(3.2, nmut), factor = 3 / 3), m3,
          col = alpha('red', 0.2), pch = 16)
  points(c(0.8, 1.8, 2.8), c(s1[nseg], s2[nseg], s3[nseg]),
          col = alpha('green', 0.75), pch = 16)
  arrows(c(0.7, 1.1, 1.7, 2.1, 2.7, 3.1),
        c(median(s1, na.rm = T), median(m1, na.rm = T), median(s2, na.rm = T), median(m2, na.rm = T), median(s3, na.rm = T), median(m3, na.rm = T)),
        c(0.9, 1.3, 1.9, 2.3, 2.9, 3.3),
        c(median(s1, na.rm = T), median(m1, na.rm = T), median(s2, na.rm = T), median(m2, na.rm = T), median(s3, na.rm = T), median(m3, na.rm = T)), length = 0)
  for (n in 1:3) {
    if (t[n] <= (0.05 / 3)) {
      text(x = n, y = 1.95, labels = signif(t[n], digits = 2), font = 4)
    } else {
      text(x = n, y = 1.95, labels = signif(t[n], digits = 2))
    }
  }
  for (i in 1:3) {
    a <- args[i]
    args[i] <- paste0('_', unlist(strsplit(a, split = '_'))[2], '_')
    if (i == 1) {
      conds <- short.ls[[which(endsWith(names(short.ls), args[i]))]]
    } else {
      conds <- c(conds, short.ls[[which(endsWith(names(short.ls), args[i]))]])
    }
  }
  comps <- c(paste0(conds[1], ':', conds[2]),
              paste0(conds[2], ':', conds[3]),
              paste0(conds[3], ':', conds[1]))
  axis(side = 1, at = c(1:3), labels = comps)
  mtext(parse(text = sprintf('italic(%s)', prom.th[pr])), side = 3, line = 0.5)
  title(ylab = 'Plasticity (in 2D)', line = 2.5)

dev.off()

### produce Extended Data Figure 4 (comparing expression values
### of all promoters between all pairs of conditions)
cat(paste('Producing Extended Data Figure 4\n'))
pdf(file = 'ED_Figure_4.pdf', width = 12, height = 12)
  par(mfrow = c(5, 6),
      las = 1,
      par(mar = c(4, 3.5, 2, 0.5)))

  for (pr in names(prom.th)) {
    args <- offnames[which(grepl(pr, offnames))]
    for (n in 1:3) {
      a <- args[n]
      if (n < 3) {
        m <- n + 1
        b <- args[m]
      } else {
        m <- 1
        b <- args[m]
      }
      mut1 <- mut.ls[[a]][, 3]
      mut1 <- as.numeric(mut1) - offset.m[a]
      mut2 <- mut.ls[[b]][, 3]
      mut2 <- as.numeric(mut2) - offset.m[b]
      seg1 <- seg.ls[[a]][, 1]
      seg2 <- seg.ls[[b]][, 1]
      nseg <- length(seg1)
      plot(mut1, mut2, pch = 16, col = alpha('red', 0.2),
          ylim = c(2, 5), xlim = c(2, 5), xlab = '', ylab = '')
      points(seg1[1:nseg - 1], seg2[1:nseg - 1],
            pch = 16, col = alpha('black', 0.2))
      points(seg1[nseg], seg2[nseg],
            pch = 16, col = alpha('green', 0.75))
      for (i in 1:3) {
        a <- args[i]
        arg <- paste0('_', unlist(strsplit(a, split = '_'))[2], '_')
        if (i == 1) {
          conds <- cond.ls[[which(endsWith(names(cond.ls), arg))]]
        } else {
          conds <- c(conds, cond.ls[[which(endsWith(names(cond.ls), arg))]])
        }
      }
      abline(0, 1, lty = 3, col = 'blue')
      mtext(parse(text = sprintf('italic(%s)', prom.th[pr])), side = 3, line = 0.5)
      title(xlab = conds[n], line = 2)
      title(ylab = conds[m], line = 2.5)
    }
  }
  legend('topleft', legend = c('Mutagenesis', 'Segregating', 'MG1655'),
          pch = 16, col = c(alpha('red', 0.2), alpha('black', 0.2), alpha('green', 0.75)),
          title = 'Promoter set')

dev.off()

##################################################
################# MEAN-CV PLOTS ##################
############# ED FIGURE 5 & FIGURE 6 #############
##################################################

### define colors for ploting and variable 'noi'
### to save noise levels in it
cols <- c('red', 'blue', 'gold')
noi <- list()

### produce Extended Data Figure 5 (fitting smooth spline
### to coefficient of variation and modal expression data)
cat(paste('Producing Extended Data Figure 5\n'))
pdf(file = 'ED_Figure_5.pdf', width = 9, height = 7)
  par(mfrow = c(3, 4),
      las = 1)

  ### loop through all promoters
  for (pr in names(prom.th)) {
    if (pr == names(prom.th)[1] || pr == names(prom.th)[5] || pr == names(prom.th)[9]) {
      par(mar = c(3.5, 5, 1.5, 0))
    } else if (pr == names(prom.th)[2] || pr == names(prom.th)[6] || pr == names(prom.th)[10]) {
      par(mar = c(3.5, 4, 1.5, 1))
    } else if (pr == names(prom.th)[3] || pr == names(prom.th)[7]) {
      par(mar = c(3.5, 3, 1.5, 2))
    } else {
      par(mar = c(3.5, 2, 1.5, 3))
    }
    ### extract info about all three environments
    args <- offnames[which(grepl(pr, offnames))]
    names(cols) <- args
    ### loop through all the environments
    for (a in args) {
      ### get and correct expression level and
      ### coefficient of variation values for random variants
      mut <- mut.ls[[a]][, 3:4]
      mut[, 1] <- as.numeric(mut[, 1]) - offset.m[a]
      mut[, 2] <- as.numeric(mut[, 2]) - offset.n[a]
      ### check for NAs in expression level and coefficient of variation
      ex1 <- as.numeric(rownames(mut)[which(is.na(mut[, 1]))])
      ex2 <- as.numeric(rownames(mut)[which(is.na(mut[, 2]))])
      ### collapse info about NAs for variants if it
      ### is present in both expression and coefficient of variation
      ex <- ex1
      if (length(ex2) > 0) {
        for (e2 in ex2) {
          hit <- 0
          for (e1 in ex) {
            if (e1 == e2) {
              hit <- 1
            }
          }
          if (hit == 0) {
            ex <- c(ex, e2)
          }
        }
      }
      ex <- sort(ex)
      ### create expression and coefficient of variation matrix
      ### without NAs for smooth spline calculation
      mut <- mut[complete.cases(mut),]
      ### calculate the smooth spline
      ### remove MG1655 variant of mtr promoter from calculation (SNP in GFP)
      if (grepl('Mtr', a)) {
        li <- smooth.spline(c(mut[, 1], seg.ls[[a]][1:(length(seg.ls[[a]][, 1]) - 1), 1]),
                            c(mut[, 2], seg.ls[[a]][1:(length(seg.ls[[a]][, 2]) - 1), 2]), lambda = 0.01)
      } else {
        li <- smooth.spline(c(mut[, 1], seg.ls[[a]][, 1]), c(mut[, 2], seg.ls[[a]][, 2]), lambda = 0.01)
      }
      ### predict coefficient of variation levels for each observed
      ### expression level from the calculated smooth spline
      fit <- predict(li, c(mut[, 1], seg.ls[[a]][, 1]))
      noise <- c(mut[, 2], seg.ls[[a]][, 2]) - fit$y
      ### add NAs into the vector of noise levels
      ### into the same places they were removed from
      ### before smooth spline calculation
      if (length(ex) > 0) {
        for (e in ex) {
          hold <- noise
          noise <- c(hold[1:(e - 1)], NA, hold[e:length(hold)])
        }
      }
      noi[[a]] <- noise
      ### plotting
      if (a == args[1]) {
        plot(mut[, 1], mut[, 2], pch = 16,
              col = alpha(cols[a], 0.3), xlim = c(2, 5), ylim = c(0.02, 0.5),
              xlab = '', ylab = '', log = 'y',
              main = parse(text = sprintf('italic(%s)', prom.th[pr])))
      } else {
        points(mut[, 1], mut[, 2], pch = 16, col = alpha(cols[a], 0.3))
      }
      if (grepl('Mtr', a)) {
        points(seg.ls[[a]][1:(length(seg.ls[[a]][, 1]) - 1), 1],
              seg.ls[[a]][1:(length(seg.ls[[a]][, 2]) - 1), 2],
              pch = 21, bg = alpha(cols[a], 0.3), col = 'black')
      } else {
        points(seg.ls[[a]][, 1], seg.ls[[a]][, 2],
              pch = 21, bg = alpha(cols[a], 0.3), col = 'black')
      }
      lines(predict(li), col = cols[a], pch = '.')
    }
    for (i in 1:3) {
      a <- args[i]
      args[i] <- paste0('_', unlist(strsplit(a, split = '_'))[2], '_')
      if (i == 1) {
        conds <- cond.ls[[which(endsWith(names(cond.ls), args[i]))]]
      } else {
        conds <- c(conds, cond.ls[[which(endsWith(names(cond.ls), args[i]))]])
      }
    }
    if (pr == names(prom.th)[7] || pr == names(prom.th)[8] || pr == names(prom.th)[9] || pr == names(prom.th)[10]) {
      title(xlab = 'Modal expression (log10, a.u.)', line = 2)
    }
    if (pr == names(prom.th)[1] || pr == names(prom.th)[5] || pr == names(prom.th)[9]) {
      title(ylab = 'mCV (stdev/mode)', line = 3)
    }
    legend('topright', legend = conds, pch = 16,
            col = c(alpha(cols[1], 0.3), alpha(cols[2], 0.3), alpha(cols[3], 0.3)),
            cex = 0.9, title = 'Condition')
  }
  plot(NULL, axes = F, ann = F, xlim = c(0, 1), ylim = c(0, 1))
  legend('topleft', legend = c('Mutagenesis', 'Segregating'),
          pch = 21, cex = 1.5, title = "Promoter set",
          col = c('red', 'black'), pt.bg = c(alpha('red', 0.3), alpha('red', 0.3)))

dev.off()

### produce Figure 6a (comparing noise
### between seg. and random variants)
cat(paste('Producing Figure 6a\n'))
pdf(file = "Figure_6a.pdf", width = 9, height = 7)
  par(mfrow = c(3, 4),
      las = 1)

  ### loop through all promoters
  for (pr in names(prom.th)) {
    if (pr == names(prom.th)[1] || pr == names(prom.th)[5] || pr == names(prom.th)[9]) {
      par(mar = c(3, 5, 2, 0))
    } else if (pr == names(prom.th)[2] || pr == names(prom.th)[6] || pr == names(prom.th)[10]) {
      par(mar = c(3, 4, 2, 1))
    } else if (pr == names(prom.th)[3] || pr == names(prom.th)[7]) {
      par(mar = c(3, 3, 2, 2))
    } else {
      par(mar = c(3, 2, 2, 3))
    }
    ### extract info about all three environments
    args <- offnames[which(grepl(pr, offnames))]
    for (i in 1:3) {
      a <- args[i]
      arg <- paste0('_', unlist(strsplit(a, split = '_'))[2], '_')
      if (i == 1) {
        conds <- short.ls[[which(endsWith(names(short.ls), arg))]]
      } else {
        conds <- c(conds, short.ls[[which(endsWith(names(short.ls), arg))]])
      }
    }
    ### loop through each environment
    for (n in 1:3) {
      a <- args[n]
      ### obtain noise values for
      ### random and seg. variants
      nmut <- length(mut.ls[[a]][, 3])
      mut <- noi[[a]][1:nmut]
      seg <- noi[[a]][nmut + 1:length(noi[[a]])]
      seg <- seg[complete.cases(seg)]
      ### remove MG1655 variant of mtr promoter from calculation (SNP in GFP)
      if (grepl('Mtr', a)) {
        ok <- seg[1:(length(seg) - 1)]
        seg <- ok
      }
      nseg <- length(seg)
      ### test the significance of difference in noise
      ### levels between seg. and random variants
      test <- wilcox.test(x = mut, y = seg, exact = F)
      ### plotting
      if (a == args[1]) {
        plot(jitter(rep(n - 0.2, nseg - 1), factor = c(3:1)[n]), seg[1:nseg - 1],
                main = parse(text = sprintf('italic(%s)', prom.th[pr])),
                xaxt = "n", xlim = c(0.5, 3.5), ylim = c(-0.025, 0.025),#ylim = c(-0.03, 0.06),#
                col = alpha('black', 0.2), pch = 16, xlab = '', ylab = '')
      } else {
        points(jitter(rep(n - 0.2, nseg - 1), factor = c(3:1)[n]), seg[1:nseg - 1],
                col = alpha('black', 0.2), pch = 16)
      }
      points(jitter(rep(n + 0.2, nmut), factor = c(3:1)[n]), mut,
              col = alpha('red', 0.2), pch = 16)
      if (grepl('Mtr', a)) {
        points(n - 0.2, seg[nseg],
                col = alpha('black', 0.2), pch = 16)
      } else {
        points(n - 0.2, seg[nseg],
                col = alpha('green', 0.75), pch = 16)
      }
      arrows(c(n - 0.3, n + 0.1), c(median(seg, na.rm = T), median(mut, na.rm = T)),
            c(n - 0.1, n + 0.3), c(median(seg, na.rm = T), median(mut, na.rm = T)), length = 0)
      if (test$p.value <= (0.05 / 3)) {
        text(x = n, y = 0.023, labels = signif(test$p.value, digits = 2), font = 4)
      } else {
        text(x = n, y = 0.023, labels = signif(test$p.value, digits = 2))
      }
    }
    axis(side = 1, at = c(1:3), labels = conds)
    if (pr == names(prom.th)[1] || pr == names(prom.th)[5] || pr == names(prom.th)[9]) {
      title(ylab = 'Noise', line = 3.5)
    }
  }
  plot(NULL, axes = F, ann = F, xlim = c(0, 1), ylim = c(0, 1))
  legend('topleft', legend = c('Mutagenesis', 'Segregating', 'MG1655'),
          pch = 16, col = c(alpha('red', 0.2), alpha('black', 0.2), alpha('green', 0.75)),
          title = 'Promoter set', cex = 1.5)

dev.off()

### produce Figure 6b (comparing noise levels
### in aceB between seg. and random variants)
cat(paste('Producing Figure 6b\n'))
pdf(file = "Figure_6b.pdf", width = 5, height = 5)
  par(las = 1)

  pr <- 'AceB'
  args <- offnames[which(grepl(pr, offnames))]
  for (i in 1:3) {
    a <- args[i]
    arg <- paste0('_', unlist(strsplit(a, split = '_'))[2], '_')
    if (i == 1) {
      conds <- short.ls[[which(endsWith(names(short.ls), arg))]]
    } else {
      conds <- c(conds, short.ls[[which(endsWith(names(short.ls), arg))]])
    }
  }
  for (n in 1:3) {
    a <- args[n]
    nmut <- length(mut.ls[[a]][, 3])
    mut <- noi[[a]][1:nmut]
    seg <- noi[[a]][nmut + 1:length(noi[[a]])]
    seg <- seg[complete.cases(seg)]
    nseg <- length(seg)
    test <- wilcox.test(x = mut, y = seg, exact = F)
    if (a == args[1]) {
      plot(jitter(rep(n - 0.2, nseg - 1), factor = c(3:1)[n]), seg[1:nseg - 1],
              main = parse(text = sprintf('italic(%s)', prom.th[pr])),
              xaxt = "n", xlim = c(0.5, 3.5), ylim = c(-0.075, 0.3),
              col = alpha('black', 0.2), pch = 16, xlab = '', ylab = '')
    } else {
      points(jitter(rep(n - 0.2, nseg - 1), factor = c(3:1)[n]), seg[1:nseg - 1],
              col = alpha('black', 0.2), pch = 16)
    }
    points(jitter(rep(n + 0.2, nmut), factor = c(3:1)[n]), mut,
            col = alpha('red', 0.2), pch = 16)
    points(n - 0.2, seg[nseg],
            col = alpha('green', 0.75), pch = 16)
    arrows(c(n - 0.3, n + 0.1), c(median(seg, na.rm = T), median(mut, na.rm = T)),
          c(n - 0.1, n + 0.3), c(median(seg, na.rm = T), median(mut, na.rm = T)), length = 0)
    if (test$p.value <= (0.05 / 3)) {
      text(x = n, y = 0.285, labels = signif(test$p.value, digits = 2), font = 4)
    } else {
      text(x = n, y = 0.285, labels = signif(test$p.value, digits = 2))
    }
  }
  axis(side = 1, at = c(1:3), labels = conds)
  title(ylab = 'Noise', line = 2.5)

dev.off()

### produce Figure 6c (comparing noise levels
### in purA between seg. and random variants)
cat(paste('Producing Figure 6c\n'))
pdf(file = "Figure_6c.pdf", width = 5, height = 5)
  par(las = 1)

  pr <- 'PurA'
  args <- offnames[which(grepl(pr, offnames))]
  for (i in 1:3) {
    a <- args[i]
    arg <- paste0('_', unlist(strsplit(a, split = '_'))[2], '_')
    if (i == 1) {
      conds <- short.ls[[which(endsWith(names(short.ls), arg))]]
    } else {
      conds <- c(conds, short.ls[[which(endsWith(names(short.ls), arg))]])
    }
  }
  for (n in 1:3) {
    a <- args[n]
    nmut <- length(mut.ls[[a]][, 3])
    mut <- noi[[a]][1:nmut]
    seg <- noi[[a]][nmut + 1:length(noi[[a]])]
    seg <- seg[complete.cases(seg)]
    nseg <- length(seg)
    test <- wilcox.test(x = mut, y = seg, exact = F)
    if (a == args[1]) {
      plot(jitter(rep(n - 0.2, nseg - 1), factor = c(3:1)[n]), seg[1:nseg - 1],
              main = parse(text = sprintf('italic(%s)', prom.th[pr])),
              xaxt = "n", xlim = c(0.5, 3.5), ylim = c(-0.025, 0.09),
              col = alpha('black', 0.2), pch = 16, xlab = '', ylab = '')
    } else {
      points(jitter(rep(n - 0.2, nseg - 1), factor = c(3:1)[n]), seg[1:nseg - 1],
              col = alpha('black', 0.2), pch = 16)
    }
    points(jitter(rep(n + 0.2, nmut), factor = c(3:1)[n]), mut,
            col = alpha('red', 0.2), pch = 16)
    points(n - 0.2, seg[nseg],
            col = alpha('green', 0.75), pch = 16)
    arrows(c(n - 0.3, n + 0.1), c(median(seg, na.rm = T), median(mut, na.rm = T)),
          c(n - 0.1, n + 0.3), c(median(seg, na.rm = T), median(mut, na.rm = T)), length = 0)
    if (test$p.value <= (0.05 / 3)) {
      text(x = n, y = 0.085, labels = signif(test$p.value, digits = 2), font = 4)
    } else {
      text(x = n, y = 0.085, labels = signif(test$p.value, digits = 2))
    }
  }
  axis(side = 1, at = c(1:3), labels = conds)
  title(ylab = 'Noise', line = 3)

dev.off()

##################################################
################### Z-SCORES #####################
################## ED FIGURE 6 ###################
##################################################

### define points in 3D on an isospline to use when
### calculating the distance of datapoints from the isospline
x1 <- rep(2, 3)
x2 <- rep(5, 3)

### produce Extended Data Figure 6 (comparing summed
### z-scores between seg. and random variants)
cat(paste('Producing Extended Data Figure 6\n'))
pdf(file = 'ED_Figure_6.pdf', width = 9, height = 5)
par(fig = c(0, 1, 0, 1),
    las = 1)

    ### loop through all promoters
    for (n in 1:length(names(prom.th))) {
      pr <- names(prom.th)[n]
      ### extract info about all three environments
      args <- offnames[which(grepl(pr, offnames))]
      a <- args[1]
      b <- args[2]
      c <- args[3]
      ### get and correct expression values
      ### from random variants
      mut1 <- mut.ls[[a]][, 3]
      mut1 <- as.numeric(mut1) - offset.m[a]
      mut2 <- mut.ls[[b]][, 3]
      mut2 <- as.numeric(mut2) - offset.m[b]
      mut3 <- mut.ls[[c]][, 3]
      mut3 <- as.numeric(mut3) - offset.m[c]
      mut <- c()
      ### calculate 3D plasticity for random variants
      for (m in 1:length(mut1)) {
        x0 <- c(mut1[m], mut2[m], mut3[m])
        mut <- c(mut, dist3d(x0, x1, x2))
      }
      ### get expression values from seg. variants
      seg1 <- seg.ls[[a]][, 1]
      seg2 <- seg.ls[[b]][, 1]
      seg3 <- seg.ls[[c]][, 1]
      ### remove MG1655 variant of mtr promoter from calculation (SNP in GFP)
      if (grepl('Mtr', pr)) {
        ok <- seg1[1:(length(seg1) - 1)]
        seg1 <- ok
        ok <- seg2[1:(length(seg2) - 1)]
        seg2 <- ok
        ok <- seg3[1:(length(seg3) - 1)]
        seg3 <- ok
      }
      seg <- c()
      ### calculate 3D plasticity for seg. variants
      for (s in 1:length(seg1)) {
        x0 <- c(seg1[s], seg2[s], seg3[s])
        seg <- c(seg, dist3d(x0, x1, x2))
      }
      ### get noise levels for both seg. and random variants
      nmut <- length(mut.ls[[a]][, 3])
      mutN1 <- noi[[a]][1:nmut]
      segN1 <- noi[[a]][nmut + 1:length(noi[[a]])]
      segN1 <- segN1[complete.cases(segN1)]
      mutN2 <- noi[[b]][1:nmut]
      segN2 <- noi[[b]][nmut + 1:length(noi[[b]])]
      segN2 <- segN2[complete.cases(segN2)]
      mutN3 <- noi[[c]][1:nmut]
      segN3 <- noi[[c]][nmut + 1:length(noi[[c]])]
      segN3 <- segN3[complete.cases(segN3)]
      ### remove MG1655 variant of mtr promoter from calculation (SNP in GFP)
      if (grepl('Mtr', pr)) {
        ok <- segN1[1:(length(segN1) - 1)]
        segN1 <- ok
        ok <- segN2[1:(length(segN2) - 1)]
        segN2 <- ok
        ok <- segN3[1:(length(segN3) - 1)]
        segN3 <- ok
      }
      ### calculate expression z-scores for 1st environment
      seg.sd <- sd(seg1)
      seg.md <- mean(seg1)
      segZ1 <- abs(seg1 - seg.md) / seg.sd
      mutZ1 <- abs(mut1 - seg.md) / seg.sd
      ### calculate expression z-scores for 2nd environment
      seg.sd <- sd(seg2)
      seg.md <- mean(seg2)
      segZ2 <- abs(seg2 - seg.md) / seg.sd
      mutZ2 <- abs(mut2 - seg.md) / seg.sd
      ### calculate expression z-scores for 3rd environment
      seg.sd <- sd(seg3)
      seg.md <- mean(seg3)
      segZ3 <- abs(seg3 - seg.md) / seg.sd
      mutZ3 <- abs(mut3 - seg.md) / seg.sd
      ### calculate plasticity z-scores
      seg.sd <- sd(seg)
      seg.md <- mean(seg)
      segZ4 <- abs(seg - seg.md) / seg.sd
      mutZ4 <- abs(mut - seg.md) / seg.sd
      ### calculate noise z-scores for 1st environment
      seg.sd <- sd(segN1)
      seg.md <- mean(segN1)
      segZ5 <- abs(segN1 - seg.md) / seg.sd
      mutZ5 <- abs(mutN1 - seg.md) / seg.sd
      ### calculate noise z-scores for 2nd environment
      seg.sd <- sd(segN2)
      seg.md <- mean(segN2)
      segZ6 <- abs(segN2 - seg.md) / seg.sd
      mutZ6 <- abs(mutN2 - seg.md) / seg.sd
      ### calculate noise z-scores for 3rd environment
      seg.sd <- sd(segN3)
      seg.md <- mean(segN3)
      segZ7 <- abs(segN3 - seg.md) / seg.sd
      mutZ7 <- abs(mutN3 - seg.md) / seg.sd
      ### check which modal expression z-scores to use
      ### (median seg. expression above 2.5 +
      ### low correlation with plasticity)
      modZS <- list('1' = segZ1, '2' = segZ2, '3' = segZ3)
      modZM <- list('1' = mutZ1, '2' = mutZ2, '3' = mutZ3)
      meds <- c(median(seg1), median(seg2), median(seg3))
      cors <- c()
      for (m in 1:3) {
        if (meds[m] > 2.5) {
          c <- cor.test(c(modZS[[m]], modZM[[m]]), c(segZ4, mutZ4), method = 'spearman', exact = F)
          cors <- c(cors, c$estimate)
        } else {
          cors <- c(cors, NA)
        }
      }
      names(cors) <- c(1:3)
      good <- names(cors)[which(!is.na(cors))]
      if (length(good) == 1) {
        use <- as.numeric(good)
      } else {
        good <- names(cors)[which(cors == min(cors, na.rm = T))]
        use <- as.numeric(good)
      }
      ### calculate summed z-scores for seg. variants
      segZ <- vector()
      for (i in 1:length(segZ1)) {
        segZ <- c(segZ, (modZS[[use]][i] + segZ4[i] + segZ5[i] + segZ6[i] + segZ7[i]))
      }
      nseg <- length(segZ)
      ### calculate summed z-scores for random variants
      mutZ <- vector()
      for (i in 1:length(mutZ1)) {
        mutZ <- c(mutZ, (modZM[[use]][i] + mutZ4[i] + mutZ5[i] + mutZ6[i] + mutZ7[i]))
      }
      nmut <- length(mutZ)
      ### test the significance of differences in summed
      ### z-scores between seg. and random variants
      test <- wilcox.test(x = mutZ, y = segZ, exact = F)
      ### plotting
      if (n == 1) {
        plot(jitter(rep(n - 0.2, nseg - 1), factor = 4 / n), segZ[1:nseg - 1],
                xaxt = 'n', ylim = c(0, 150), xlim = c(0.5, 10.5),
                col = alpha('black', 0.2), pch = 16, xlab = '', ylab = '')
      } else {
        points(jitter(rep(n - 0.2, nseg - 1), factor = 4 / n), segZ[1:nseg - 1],
                col = alpha('black', 0.2), pch = 16)
      }
      points(jitter(rep(n + 0.2, nmut), factor = 4 / n), mutZ,
              col = alpha('red', 0.2), pch = 16)
      if (grepl('Mtr', pr)) {
        points(n - 0.2, segZ[nseg],
                col = alpha('black', 0.2), pch = 16)
      } else {
        points(n - 0.2, segZ[nseg],
                col = alpha('green', 0.75), pch = 16)
      }
      arrows(c(n - 0.3, n + 0.1), c(median(segZ, na.rm = T), median(mutZ, na.rm = T)),
            c(n - 0.1, n + 0.3), c(median(segZ, na.rm = T), median(mutZ, na.rm = T)), length = 0)
      if (test$p.value < 0.05) {
        text(x = n, y = 150, labels = signif(test$p.value, digits = 2), font = 4)
      } else {
        text(x = n, y = 150, labels = signif(test$p.value, digits = 2))
      }
    }
    axis(side = 1, at = c(1:10), labels = parse(text = sprintf('italic(%s)', prom.th)))
    title(ylab = 'Cumulative z-score', line = 2.5)
    legend(0.1, 70, legend = c('Mutagenesis', 'Segregating', 'MG1655'),
            pch = 16, col = c(alpha('red', 0.2), alpha('black', 0.2), alpha('green', 0.75)),
            title = 'Promoter set', cex = 0.75)

### plotting inset
par(fig = c(0.05, 0.35, 0.37, 0.94),
    new = TRUE)
  for (i in 1:length(names(prom.th))) {
    pr <- names(prom.th)[i]
    if (grepl('AceB', pr) || grepl('PurA', pr)) {
      args <- offnames[which(grepl(pr, offnames))]
      a <- args[1]
      b <- args[2]
      c <- args[3]
      mut1 <- mut.ls[[a]][, 3]
      mut1 <- as.numeric(mut1) - offset.m[a]
      mut2 <- mut.ls[[b]][, 3]
      mut2 <- as.numeric(mut2) - offset.m[b]
      mut3 <- mut.ls[[c]][, 3]
      mut3 <- as.numeric(mut3) - offset.m[c]
      mut <- c()
      x1 <- rep(2, 3)
      x2 <- rep(5, 3)
      for (m in 1:length(mut1)) {
        x0 <- c(mut1[m], mut2[m], mut3[m])
        mut <- c(mut, dist3d(x0, x1, x2))
      }

      seg1 <- seg.ls[[a]][, 1]
      seg2 <- seg.ls[[b]][, 1]
      seg3 <- seg.ls[[c]][, 1]
      seg <- c()
      for (s in 1:length(seg1)) {
        x0 <- c(seg1[s], seg2[s], seg3[s])
        seg <- c(seg, dist3d(x0, x1, x2))
      }

      nmut <- length(mut.ls[[a]][, 3][which(complete.cases(mut.ls[[a]]))])
      mutN1 <- noi[[a]][1:nmut]
      segN1 <- noi[[a]][nmut + 1:length(noi[[a]])]
      segN1 <- segN1[complete.cases(segN1)]
      nmut <- length(mut.ls[[b]][, 3][which(complete.cases(mut.ls[[b]]))])
      mutN2 <- noi[[b]][1:nmut]
      segN2 <- noi[[b]][nmut + 1:length(noi[[b]])]
      segN2 <- segN2[complete.cases(segN2)]
      nmut <- length(mut.ls[[c]][, 3][which(complete.cases(mut.ls[[c]]))])
      mutN3 <- noi[[c]][1:nmut]
      segN3 <- noi[[c]][nmut + 1:length(noi[[c]])]
      segN3 <- segN3[complete.cases(segN3)]

      seg.sd <- sd(seg1)
      seg.md <- mean(seg1)
      segZ1 <- abs(seg1 - seg.md) / seg.sd
      mutZ1 <- abs(mut1 - seg.md) / seg.sd

      seg.sd <- sd(seg2)
      seg.md <- mean(seg2)
      segZ2 <- abs(seg2 - seg.md) / seg.sd
      mutZ2 <- abs(mut2 - seg.md) / seg.sd

      seg.sd <- sd(seg3)
      seg.md <- mean(seg3)
      segZ3 <- abs(seg3 - seg.md) / seg.sd
      mutZ3 <- abs(mut3 - seg.md) / seg.sd

      seg.sd <- sd(seg)
      seg.md <- mean(seg)
      segZ4 <- abs(seg - seg.md) / seg.sd
      mutZ4 <- abs(mut - seg.md) / seg.sd

      seg.sd <- sd(segN1)
      seg.md <- mean(segN1)
      segZ5 <- abs(segN1 - seg.md) / seg.sd
      mutZ5 <- abs(mutN1 - seg.md) / seg.sd

      seg.sd <- sd(segN2)
      seg.md <- mean(segN2)
      segZ6 <- abs(segN2 - seg.md) / seg.sd
      mutZ6 <- abs(mutN2 - seg.md) / seg.sd

      seg.sd <- sd(segN3)
      seg.md <- mean(segN3)
      segZ7 <- abs(segN3 - seg.md) / seg.sd
      mutZ7 <- abs(mutN3 - seg.md) / seg.sd

      segZ <- vector()
      for (i in 1:length(segZ1)) {
        segZ <- c(segZ, (segZ1[i] + segZ2[i] + segZ3[i] + segZ4[i] + segZ5[i] + segZ6[i] + segZ7[i]))
      }
      nseg <- length(segZ)

      mutZ <- vector()
      for (i in 1:length(mutZ1)) {
        mutZ <- c(mutZ, (mutZ1[i] + mutZ2[i] + mutZ3[i] + mutZ4[i] + mutZ5[i] + mutZ6[i] + mutZ7[i]))
      }
      nmut <- length(mutZ)

      if (grepl('AceB', pr)) {
        n <- 1
        plot(jitter(rep(n - 0.2, nseg - 1), factor = 4 / n), segZ[1:nseg - 1],
                xaxt = 'n', ylim = c(0, 390), xlim = c(0.5, 2.5), cex.axis = 0.7,
                col = alpha('black', 0.2), pch = 16, xlab = '', ylab = '', cex = 0.5)
      } else {
        n <- 2
        points(jitter(rep(n - 0.2, nseg - 1), factor = 4 / n), segZ[1:nseg - 1],
                col = alpha('black', 0.2), pch = 16, cex = 0.5)
      }
      points(jitter(rep(n + 0.2, nmut), factor = 4 / n), mutZ,
              col = alpha('red', 0.2), pch = 16, cex = 0.5)
      points(n - 0.2, segZ[nseg],
              col = alpha('green', 0.75), pch = 16, cex = 0.5)
      arrows(c(n - 0.3, n + 0.1), c(median(segZ, na.rm = T), median(mutZ, na.rm = T)),
            c(n - 0.1, n + 0.3), c(median(segZ, na.rm = T), median(mutZ, na.rm = T)), length = 0)
    }
  }
  axis(side = 1, padj = -1, at = c(1:2), labels = parse(text = sprintf('italic(%s)', c('aceB', 'purA'))), cex.axis = 0.7)

dev.off()

##################################################
################## CELL GATING ###################
################## ED FIGURE 7 ###################
##################################################

### produce Extended Data Figure 7 (gating example)
cat(paste('Producing Extended Data Figure 7\n'))
pdf(file = "ED_Figure_7.pdf", width = 12, height = 3)
  par(las = 1, xpd = NA,
      mar = c(4, 4, 2, 1),
      mfrow = c(1, 4))

    ### get all files names ending with .fsc
    Files <- Sys.glob(paste0(data.path[1], '*.fcs'))
    ### read sample data file
    data <- read.FCS(Files[5], alter.names = T)
    ### indentify values for xlim and ylim
    xs <- c(min(log10(as.vector(exprs(data$FSC.H)))),
            max(log10(as.vector(exprs(data$FSC.H)))))
    ys <- c(min(log10(as.vector(exprs(data$SSC.H)))),
            max(log10(as.vector(exprs(data$SSC.H)))))
    ### plot raw data of forward and side scatter
    plot(log10(as.vector(exprs(data$FSC.H))),
          log10(as.vector(exprs(data$SSC.H))),
          main = 'Raw Data', xlab = '', ylab = '',
          pch = '.', xlim = xs, ylim = ys,
          col = alpha('black', 0.1))
    title(xlab = 'FSC.H (log10, a.u.)', line = 2.5)
    title(ylab = 'SSC.H (log10, a.u.)', line = 2.5)
    text(2, 5.6, 'a', cex = 1.5, font = 2)
    mtext(text = paste(' n =', length(as.vector(exprs(data$FSC.H)))),
          side = 3, adj = 0, line = -1, cex = 0.75)

    ### get max kernel density values from scatter values
    fsc.dens <- density(as.vector(exprs(data$FSC.H)))
    ssc.dens <- density(as.vector(exprs(data$SSC.H)))
    fypeak <- which.max(fsc.dens$y)
    sypeak <- which.max(ssc.dens$y)
    f.peak <- fsc.dens$x[fypeak]
    s.peak <- ssc.dens$x[sypeak]
    ### plot raw data with the highest kernel density point
    plot(log10(as.vector(exprs(data$FSC.H))),
          log10(as.vector(exprs(data$SSC.H))),
          main = 'Highest kernel density detection', xlab = '', ylab = '',
          pch = '.', xlim = xs, ylim = ys,
          col = alpha('black', 0.1))
    title(xlab = 'FSC.H (log10, a.u.)', line = 2.5)
    title(ylab = 'SSC.H (log10, a.u.)', line = 2.5)
    text(2, 5.6, 'b', cex = 1.5, font = 2)
    mtext(text = paste(' n =', length(as.vector(exprs(data$FSC.H)))),
          side = 3, adj = 0, line = -1, cex = 0.75)
    points(log10(f.peak), log10(s.peak), pch = '+',
            col = 'red')

    ### create rectangle gate at based on the kernel density
    pregate <- rectangleGate("FSC.H" = c(f.peak - 2000, f.peak + 15000),
                            "SSC.H" = c(s.peak - 2000, s.peak + 6000))
    filt <- filter(data, pregate)
    data <- Subset(data, filt)
    ### plot remaining events
    plot(log10(as.vector(exprs(data$FSC.H))),
          log10(as.vector(exprs(data$SSC.H))),
          main = 'Outlier removal', xlab = '', ylab = '',
          pch = '.', xlim = xs, ylim = ys,
          col = alpha('black', 0.1))
    title(xlab = 'FSC.H (log10, a.u.)', line = 2.5)
    title(ylab = 'SSC.H (log10, a.u.)', line = 2.5)
    text(2, 5.6, 'c', cex = 1.5, font = 2)
    mtext(text = paste(' n =', length(as.vector(exprs(data$FSC.H)))),
          side = 3, adj = 0, line = -1, cex = 0.75)
    points(log10(f.peak), log10(s.peak), pch = '+',
            col = 'red')

    ### then create an elipse gate from cells filtered by rectangle gate
    FSC.H <- as.vector(exprs(data$FSC.H))
    SSC.H <- as.vector(exprs(data$SSC.H))
    SC.H <- cbind(FSC.H, SSC.H)
    cv <- cov(SC.H)                         ### covariance matrix needed for elipse gate definition
    mn <- c(median(FSC.H), median(SSC.H))   ### mean to define the elipse gate
    elgate <- ellipsoidGate(.gate = cv, mean = mn, distance = 1)
    filt <- filter(data, elgate)
    cells <- Subset(data, filt)
    ### plot events passing through the final gate
    plot(log10(as.vector(exprs(cells$FSC.H))),
          log10(as.vector(exprs(cells$SSC.H))),
          main = 'Core density gating', xlab = '', ylab = '',
          pch = '.', xlim = xs, ylim = ys,
          col = alpha('black', 0.1))
    title(xlab = 'FSC.H (log10, a.u.)', line = 2.5)
    title(ylab = 'SSC.H (log10, a.u.)', line = 2.5)
    text(2, 5.6, 'd', cex = 1.5, font = 2)
    mtext(text = paste(' n =', length(as.vector(exprs(cells$FSC.H)))),
          side = 3, adj = 0, line = -1, cex = 0.75)
    points(log10(f.peak), log10(s.peak), pch = '+',
            col = 'red')

dev.off()
