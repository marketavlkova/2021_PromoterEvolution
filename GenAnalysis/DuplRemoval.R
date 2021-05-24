#!/usr/bin/env Rscript
### Script for removing duplicated sequences in a fasta file
### example of use: ./DuplRemoval.R <input file> <output file>

### read arguments from command line
args = commandArgs(trailingOnly = TRUE)

### load data and calculate its length
data <- read.csv(args[1], header = F)
len <- length(data[, 1])

### set output variable 'clear' for sequences without duplication
### loop through all sequences
clear <- 1
for (j in seq(1, len - 1, by = 2)) {
  ### if some promoter is already saved in 'clear' variable
  if (!is.numeric(clear)) {
    name <- as.character(data[j, 1])
    prom <- unlist(strsplit(name, split = '_'))[4]
    ### if previous and current promoter names do not match
    ### add it into the 'clear' variable
    if (!grepl(b, prom)) {
      clear <- rbind(clear, as.character(data[j, 1]))
      clear <- rbind(clear, as.character(data[j + 1, 1]))
      a <- unlist(strsplit(prom, ''))
      join <- a[1]
      for (i in 2:length(a)) {
        if (!grepl('[0-9]', a[i]) && !grepl('-', a[i])) {
          join <- paste0(join, a[i])
        } else {
          break
        }
      }
      b <- join
      names <- c(names, b)
    }
  ### in case of first round through the loop
  ### overwrite the 'clear' variable with the current sequence
  ### and save its promoter name to check against in the nest round
  } else {
    name <- as.character(data[j, 1])
    prom <- unlist(strsplit(name, split = '_'))[4]
    a <- unlist(strsplit(prom, ''))
    join <- a[1]
    for (i in 2:length(a)) {
      if (!grepl('[0-9]', a[i]) && !grepl('-', a[i])) {
        join <- paste0(join, a[i])
      } else {
        break
      }
    }
    b <- join
    names <- b
    clear <- rbind(as.character(data[j, 1]), as.character(data[j + 1, 1]))
  }
}

### generate output file which does not contran duplicates
write.table(clear, file = args[2], sep = ',', col.names = FALSE, row.names = FALSE, quote = FALSE)
