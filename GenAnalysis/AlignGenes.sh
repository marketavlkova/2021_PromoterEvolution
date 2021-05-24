#!/bin/bash
### Script for alignment and average pairwise identity calculation of genes among strains
### example to run: bash AlignGenes.sh <path to directory with fasta files containing sequences to be aligned>

### loop through all fasta files and perform multiple sequence alignment
for file in $1*.fasta
do
  NO_SUFIX="${file%.*}"
  NAME="${NO_SUFIX#*/*/*/}"
  echo "Aligning $NAME"
  t_coffee $file -type=dna -run_name=output/genes/$NAME &> /dev/null
done

### delete file with average pairwise identities if it already exists
if [ -f "output/GeneIdentities.txt" ]
then
  rm output/GeneIdentities.txt
fi

### calculate average pairwise identities for aligned sequences
echo "Calculating identities for aligned sequences"
for align in output/genes/*.aln
do
  t_coffee -other_pg seq_reformat -in $align -output sim 2> /dev/null
done >> output/GeneIdentities.txt
