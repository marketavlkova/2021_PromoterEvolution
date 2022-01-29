#!/bin/bash
### Workflow script for sequence analysis
### optimized to run on MacOS High Sierra 10.13.6, using: Awk 20070501, R 4.0.3
### created by Marketa Vlkova in 24-May-2021
### example to run: ./RunAnalysis.sh

### STEP 1:
### get files of aligned promoters having each sequence on one line only in txt files
printf "Producing text files from alignments for variant number calculation\n"
for file in output/promoters/*.aln
do
  OUT="${file%.aln}.txt"
  ./JoinLines.awk $file > $OUT
done
for file in output/IGRs/*.aln
do
  OUT="${file%.aln}.txt"
  ./JoinLines.awk $file > $OUT
done
for file in output/genes/*.aln
do
  OUT="${file%.aln}.txt"
  ./JoinLines.awk $file > $OUT
done
### get number of sequence variants for alignment
printf "Pulling out numbers of existing versions for:\n"
printf "\tpromoters\n"
cat output/promoters/*.txt | ./PromVerCount.awk > "output/VariantProms.csv"
printf "\tIGRs\n"
cat output/IGRs/*.txt | ./PromVerCount.awk > "output/VariantIGRs.csv"
printf "\tgenes\n"
cat output/genes/*.txt | ./PromVerCount.awk > "output/VariantGenes.csv"

### STEP 2:
### extract function groups from MultiFun tables
printf "Extracting data from MultiFun database files\n"
./MultiFunKey.sh 'input/MultiFun/IDs?.txt' output/MultiFunKey.csv output/MultiFunInput.csv

### STEP 3:
### calculate Pi & Watterson estimator Theta
printf "Calculating Pi and Watterson estimator Theta"
./ThetaPi.R

### STEP 4:
### produce Supp Figures (1c and 1d) checking variability
### in promoters with and without flanking ORFs
printf "Producing Supp Figures 1a and 1b\n"
./CorrTestPlots.R

### STEP 5:
### create list of gene names for function group plots
printf "Generating list of gene names for function group analysis\n"
cat output/DataIGRs\&GenesCommon.csv | awk 'BEGIN{
  FS = ",";
  printf "Gene";
}
{
  N = split($1, L, "");
  if (L[N] ~ /^[0-9]+$/) {
    M = N - 1;
  } else if (L[N-1] ~ /^[0-9]+$/) {
    M = N - 2;
  } else if ($1 ~ /^udpP/) {
    M = N + 1;
  } else {
    M = N;
  }
  for (I = 1; I < M; I++) {
    printf "%s", L[I];
  }
  printf "\n";
}' > output/GeneNames.csv

### STEP 6:
### produce Figures (1a and Supp 1a) comparing variability
### in promoters and their downstream ORFs
printf "Producing Figure 1a and Figure 1c\n"
python3 BokehPlot.py

### STEP 7:
### produce Figures (1b and Supp 1b) checking whether
### any functional groups are enriched for promoters
### with high or low variability
printf "Producing Figure 1b and Figure 1d\n"
./FunctionPlots.R
