#!/bin/bash
### Script for alignment and average pairwise identity calculation of intergenic / promoter regions among strains
### example to run: bash Align.sh <path to directory with fasta files containing sequences to be aligned> <output directory> <number of which alignment length can be smaller than MG1655 sequence to be included>

### loop through all fasta files and perform multiple sequence alignment (promoter specific)
for file in $1*$3*.fasta
do
  NO_SUFIX="${file%.*}"
  NAME="${NO_SUFIX#*/*/*/}"
  echo "Aligning $NAME"
  t_coffee $file -mode procoffee -run_name=$2$NAME &> /dev/null
done

### delete file with average pairwise identities if it already exists
if [ -f "output/PromIdentities.txt" ]
then
  rm output/PromIdentities.txt
fi

### calculate average pairwise identities for aligned sequences
echo "Calculating identities for aligned promoters"
for align in $2*$3*.aln
do
  t_coffee -other_pg seq_reformat -in $align -output sim 2> /dev/null
done >> output/IdentitiesProms.txt
