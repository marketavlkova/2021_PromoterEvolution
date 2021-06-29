#!/bin/bash
### Workflow script for analysis of fluorescence data
### optimized to run on MacOS High Sierra 10.13.6, using: Python 3.9.1, R 4.0.3
### created by Marketa Vlkova in 29-June-2021
### example to run: ./Run.sh

### STEP 1:
### get positions of SNPs from random variants with just single SNP
printf "Producing SNP position information for random variants with single SNPs.\n"
python3 SNPmap.py

### STEP 2:
### produce Figures (3 - 8 and Supp 2 - 5)
./PlotPhenotype.R
