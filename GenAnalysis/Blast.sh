#!/bin/bash
### Script for blasting sequences obtained from MG1655 reference against other genomes
### example to run: ./Blast.sh [should be database created? Yes/No] <path to fasta file with MG1655 reference sequences> <path to output files>

### set path to blast database and output files
DB_PATH="input/Database"
OUT_PATH=$3
### set path to query file
QUERY_FILE=$2

### check whether genome database should be made
if [ "$1" != "${1#[Yy]}" ]
then
  ### loop through all fasta (.fsa) genomes provided and make database for each of them
  for FASTA in input/genomes/*/*.fsa
  do
    PRE_DB="${FASTA#*/*/*/}"
    FIN_DB="${PRE_DB%.*}"
    makeblastdb -in $FASTA -dbtype nucl -out $DB_PATH/$FIN_DB
  done
fi

### do the blasting for each genome
for FILE in $DB_PATH/*.nhr
do
  IN_DB="${FILE#*/*/}"
  OUT_FILE="${IN_DB%.*}_blasted.txt"
  printf "Blasting $IN_DB \n"
  blastn -db $DB_PATH/"${IN_DB%.*}" -query $QUERY_FILE -strand both -task megablast -evalue 1e-10 -out $OUT_PATH$OUT_FILE -outfmt "6 std" -num_threads 4
done
