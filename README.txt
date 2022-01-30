==== README for supplementary data and scripts ====

The workflow is divided into 3 main parts defined by scripts:
1. RunExtraction.sh (extracting sequences from strains, located in GenAnalysis directory)
2. RunAnalysis.sh (analysing and plotting sequence information results, located in GenAnalysis directory)
3. Run.sh (analysing and plotting phenotypic data, located in the main directory)
(run in this order)

Each of the three scripts mentioned above combines smaller parts (scripts) of the workflow in an appropriate order.
RunExtraction.sh scripts are optimised to run on Ubuntu OS (see details and dependencies in the header of the script).
RunAnalysis.sh and Run.sh scripts are optimised to run on macOS (see details and dependencies in the header of the two scripts).

Each individual script has a description of function at the beginning of the file as well as in one of the 3 workflow scripts at the an appropriate position.
All output figures will appear in the main directory, except for an interactive Bokeh plot outputs (Fig. 1a and Supplementary Figure 1a) which will appear in a browser tab, from where they can be downloaded.

===================================================

Necessary libraries and packages:
(this information can be also found at the beginning of each individual R or python script)

For R scripts:
'flowCore' package from Bioconductor for flow cytometry data analysis
'scales' library
'scatterplot3d' library

For python scripts:
'Biopython'
'bokeh'
'requests'
'pandas'
'statistics'
'matplotlib'
'io'
'sys'
'os'
'glob'
'csv'
'tqdm'
