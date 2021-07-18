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
### remove duplicated promoters and genes in reference
printf "Removing duplications\n"
./DuplRemoval.R output/*ExtProms.fasta output/Proms_refs.fasta
./DuplRemoval.R output/*genes.fasta output/Genes_refs.fasta

### STEP 4:
### get genome files of environmental strains from Figshare
cd input
wget --content-disposition https://ndownloader.figshare.com/files/28871856
tar -xzvf genomes.tar.gz
rm genomes.tar.gz
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
if [ ! -d "./output/blast_results_promoters" ]
then
  mkdir output/blast_results_promoters
fi
if [ ! -d "./output/blast_results_genes" ]
then
  mkdir output/blast_results_genes
fi
### blast for homologous sequences among environmental isolates
./Blast.sh $MAKEDB ./output/Proms_refs.fasta output/blast_results_promoters/
./Blast.sh No ./output/Genes_refs.fasta output/blast_results_genes/

### STEP 5:
### create output directories for sequences found in environmental isolates
if [ ! -d "./output/blasted_promoters" ]
then
  mkdir output/blasted_promoters
  mkdir output/blasted_promoters/by_promoter
fi
if [ ! -d "./output/blasted_genes" ]
then
  mkdir output/blasted_genes
  mkdir output/blasted_genes/by_promoter
fi
### extract blasted sequences from environmental strains
printf "\nAnalysing blast promoter hits that are at most 100nt shorter than MG1655 reference"
python3 GetBlasted.py output/Proms_refs.fasta output/blast_results_promoters/ input/genomes/ output/blasted_promoters/ 100
printf "\nAnalysing blast gene hits that are at most 100nt shorter than MG1655 reference"
python3 GetBlasted.py output/Genes_refs.fasta output/blast_results_genes/ input/genomes/ output/blasted_genes/ 100

### STEP 6:
### create output directory for alignments
if [ ! -d "./output/promoters" ]
then
  mkdir output/promoters
fi
if [ ! -d "./output/genes" ]
then
  mkdir output/genes
fi
### align homologous sequences and calculate sequence identities
printf "\nAligning:\n"
printf "\tpromoter hits\n"
./AlignProms.sh output/blasted_promoters/by_promoter/ output/promoters/ 100
printf "\tgene hits\n"
./AlignGenes.sh output/blasted_genes/by_promoter/
### extract total identity values for each alignment in a separated file (average pairwise identity)
printf "Extracting total identity values from:\n"
printf "\tpromoters\n"
./APIextract.awk output/IdentitiesProms.txt > output/APIprom.csv
printf "\tgenes\n"
./APIextract.awk output/IdentitiesGenes.txt > output/APIgene.csv

### STEP 7:
### create output directory for IGR alignments
if [ ! -d "./output/IGRs" ]
then
  mkdir output/IGRs
fi
### generate IGR alignments from promoter alignments (IGR + flanks)
python3 IGRextract.py output/promoters/ output/IGRs/
### extract total identity values for IGR alignments (average pairwise identity)
echo "Calculating identities for aligned IGRs"
for align in output/IGRs/*.aln
do
  t_coffee -other_pg seq_reformat -in $align -output sim 2> /dev/null
done >> output/IdentitiesIGRs.txt
printf "Extracting total identity values from:\n"
printf "\tIGRs\n"
./APIextract.awk output/IdentitiesIGRs.txt > output/APIigr.csv

### STEP 8:
### calculate number of segregating sites for each sequence alignment
printf "Pulling out information about segregating sites from:\n"
printf "\tpromoters\n"
python3 SegSites.py output/promoters/ output/SegSitesProms.csv
printf "\tIGRs\n"
python3 SegSites.py output/IGRs/ output/SegSitesIGRs.csv
printf "\tgenes\n"
python3 SegSites.py output/genes/ output/SegSitesGenes.csv
