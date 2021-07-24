#!/bin/bash
### Workflow script for analysis of fluorescence data
### optimized to run on MacOS High Sierra 10.13.6, using: Python 3.9.1, R 4.0.3
### created by Marketa Vlkova in 25-July-2021
### example to run: ./Run.sh

### STEP 1:
### get input files for phenotypic analysis data analysis from Figshare
mkdir PhenAnalysis
cd PhenAnalysis
wget --content-disposition https://ndownloader.figshare.com/files/28871868
wget --content-disposition https://ndownloader.figshare.com/files/28871871
wget --content-disposition https://ndownloader.figshare.com/files/28871877
wget --content-disposition https://ndownloader.figshare.com/files/28871880
wget --content-disposition https://ndownloader.figshare.com/files/28871886
wget --content-disposition https://ndownloader.figshare.com/files/28871889
wget --content-disposition https://ndownloader.figshare.com/files/28871907
wget --content-disposition https://ndownloader.figshare.com/files/28871913
wget --content-disposition https://ndownloader.figshare.com/files/28871922
wget --content-disposition https://ndownloader.figshare.com/files/28871928
tar -xzvf AceB.tar.gz
rm AceB.tar.gz
tar -xzvf AldA.tar.gz
rm AldA.tar.gz
tar -xzvf Cdd.tar.gz
rm Cdd.tar.gz
tar -xzvf DctA.tar.gz
rm DctA.tar.gz
tar -xzvf LacZ.tar.gz
rm LacZ.tar.gz
tar -xzvf Mtr.tar.gz
rm Mtr.tar.gz
tar -xzvf PtsG.tar.gz
rm PtsG.tar.gz
tar -xzvf PurA.tar.gz
rm PurA.tar.gz
tar -xzvf TpiA.tar.gz
rm TpiA.tar.gz
tar -xzvf YhjX.tar.gz
rm YhjX.tar.gz
cd -

### STEP 2:
### get positions of SNPs from random variants with just single SNP
printf "Producing SNP position information for random variants with single SNPs.\n"
python3 SNPmap.py

### STEP 3:
### produce Figures (3 - 8 and Supp 2 - 5)
./PlotPhenotype.R
