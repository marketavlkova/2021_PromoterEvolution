### Script which counts number of segregation sites for each intergenic / promoter region
### Note: This script isn't optimized to work on MacOS
### example to run: python3 IGRexctract.py <path to dir with .aln files> <path to output files>

import sys
import glob
from Bio import AlignIO
from Bio.Align import AlignInfo
from tqdm import tqdm

def main():

    ### save arguments
    in_dir = sys.argv[1]
    out_dir = sys.argv[2]

    ### read all files with blast results
    in_files = glob.glob(in_dir + '*.aln')

    print("\nExtracting IGRs from promoters:")
    pbar = tqdm(in_files)
    ### loop through all input files
    for file in pbar:
        ### load sequence alignment
        align = AlignIO.read(open(file), "clustal")
        promoter = file.split("/")[2].split("_")[0]
        pbar.set_description("Processing %s" % promoter)
        len = align.get_alignment_length()
        out_file = open(out_dir + '{}.aln'.format(promoter), 'w')
        ### generate output file without 100bp flanking regions
        if ((len - 202) > 5):
            igr = align[:, 100:(len - 102)]
            out_file.write(format(igr, "clustal"))
        ### correction for one shorter promoter
        else:
            igr = align[:, 100:(len - 94)]
            out_file.write(format(igr, "clustal"))

if __name__ == '__main__':
    main()
