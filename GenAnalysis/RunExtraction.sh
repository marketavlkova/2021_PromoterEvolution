#!/bin/bash
### Workflow script for sequence search and extraction
### optimized to run on Ubuntu 16.04.7, using: Awk 4.1.3, Python 3.6.7, R 4.0.2, BLAST 2.2.31+, T-COFFEE 11.0.8
### created by Marketa Vlkova in 24-May-2021
### example to run: ./RunExtraction.sh

### Explanations for shortcuts or other terms:
### TSS = transcription start site
### IGR = intergenic region
### flanks = IGRs + flanking regions (100 bp upstream and downstream)

### create output directory
if [ ! -d "./output" ]
then
 mkdir output
fi

### STEP 1:
### extract TSS information from desired columns in RegulonDB file
printf "Extracting TSS information from Regulon DB file\n"
./RegDBextract.awk -v COLS=2,3,4,6,8 input/MG1655_s70_promoters_RegDB.txt > output/RDBex.csv

### STEP 2:
### search for and extract reference sequences from MG1655 genome using output file from STEP 1
### use only promoters with Strong confidence level
python3 GetPromRefs.py output/RDBex.csv input/MG1655_genome.gb Strong
python3 GetGeneRefs.py output/RDBex.csv input/MG1655_genome.gb Strong

### STEP 3:
### remove duplicated sequences in reference (due to multiple TSSs)
printf "Removing duplications\n"
./DuplRemoval.R output/*ExtFlanks.fasta output/ExtFlanks_refs.fasta
./DuplRemoval.R output/*ExtIGRs.fasta output/ExtIGRs_refs.fasta
./DuplRemoval.R output/*genes.fasta output/Genes_refs.fasta

### STEP 4:
### extract zip file with genomes of env. strains
cd input
tar -xzvf genomes.tar.gz
cd -
### create input/Database directory if it doesn't exist
if [ ! -d "./input/Database" ]
then
  mkdir input/Database
  MAKEDB="Yes"
else
  MAKEDB="No"
fi
### create output directory for blast results
if [ ! -d "./output/blast_results_IGRs" ]
then
  mkdir output/blast_results_IGRs
fi
if [ ! -d "./output/blast_results_flanks" ]
then
  mkdir output/blast_results_flanks
fi
if [ ! -d "./output/blast_results_genes" ]
then
  mkdir output/blast_results_genes
fi
### blast for homologous sequences among environmental isolates
./Blast.sh $MAKEDB ./output/ExtFlanks_refs.fasta output/blast_results_flanks/
./Blast.sh No ./output/ExtIGRs_refs.fasta output/blast_results_IGRs/
./Blast.sh No ./output/Genes_refs.fasta output/blast_results_genes/

### STEP 5:
### create output directories for sequences found in environmental isolates
if [ ! -d "./output/blasted_IGRs" ]
then
  mkdir output/blasted_IGRs
  mkdir output/blasted_IGRs/by_promoter
fi
if [ ! -d "./output/blasted_flanks" ]
then
  mkdir output/blasted_flanks
  mkdir output/blasted_flanks/by_promoter
fi
if [ ! -d "./output/blasted_genes" ]
then
  mkdir output/blasted_genes
  mkdir output/blasted_genes/by_promoter
fi
### extract blasted sequences from environmental strains
printf "\nAnalysing blast flank hits that are at most 100nt shorter than MG1655 reference"
python3 GetBlasted.py output/ExtFlanks_refs.fasta output/blast_results_flanks/ input/genomes/ output/blasted_flanks/ 100
printf "\nAnalysing blast IGR hits that are at most 5nt shorter than MG1655 reference"
python3 GetBlasted.py output/ExtIGRs_refs.fasta output/blast_results_IGRs/ input/genomes/ output/blasted_IGRs/ 5
printf "\nAnalysing blast gene hits that are at most 100nt shorter than MG1655 reference"
python3 GetBlasted.py output/Genes_refs.fasta output/blast_results_genes/ input/genomes/ output/blasted_genes/ 100

### STEP 6:
### create output directory for alignments
if [ ! -d "./output/flanks" ]
then
  mkdir output/flanks
fi
if [ ! -d "./output/IGRs" ]
then
  mkdir output/IGRs
fi
if [ ! -d "./output/genes" ]
then
  mkdir output/genes
fi
### align homologous sequences and calculate sequence identities
printf "\nAligning:\n"
printf "\tflank hits\n"
./AlignProms.sh output/blasted_flanks/by_promoter/ output/flanks/ 100
printf "\tIGR hits\n"
./AlignProms.sh output/blasted_IGRs/by_promoter/ output/IGRs/ 5
printf "\tgene hits\n"
./AlignGenes.sh output/blasted_genes/by_promoter/
### extract total identity values for each alignment in a separated file (average pairwise identity)
printf "Extracting total identity values from:\n"
printf "flanks\n"
./APIextract.awk output/Identities-100.txt > output/APIprom-100.csv
printf "IGRs\n"
./APIextract.awk output/Identities-5.txt > output/APIprom-5.csv
printf "genes\n"
./APIextract.awk output/GeneIdentities.txt > output/APIgene.csv
### get files of aligned promoters having each sequence on one line only in txt files
### (these txt files are used to calculate number of sequence versions among all strains later)
for file in output/flanks/*.aln
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

### STEP 7:
### calculate number of segregating sites for each sequence alignment
printf "Pulling out information about segregating sites from:\n"
printf "flanks\n"
python3 SegSitesProms.py output/flanks/ 100
printf "IGRs\n"
python3 SegSitesProms.py output/IGRs/ 5
printf "genes\n"
python3 SegSitesGenes.py output/genes/
