### Script which counts number of segregation sites for each promoter
### example to run: python3 SNPmap.py

import sys
import glob
from Bio import AlignIO
from Bio.Align import AlignInfo

def main():

    ### save all input directories
    in_dirs = glob.glob('PhenAnalysis/*')

    ### loop through all promoter types
    for dir in in_dirs:
        ### initiate output file
        out = open(dir + '/1SNPmap.csv', 'w')
        out.write('ID,SNP_position\n')
        ### read input file
        align = AlignIO.read(open(dir + '/MutagenesisAlignment.phy'), "phylip")
        ### loop through all sequences
        for record in align:
            ### if handling MG1655 variant, save it as reference
            if ('MG1655' in record.id):
                ref = record.seq
            ### for all random variants count number of SNPs
            ### as compared to reference (MG1655 variant)
            ### and remember position of the last SNP
            else:
                snp = 0
                for i in range(0, len(ref)):
                    if (record.seq[i] != ref[i]):
                        snp += 1
                        hit = i + 1
                ### for variants with a single SNP save
                ### variant ID and SNP position in output file
                if (snp == 1):
                    out.write('%s,%s\n' % (record.id, hit))

        out.write('max,%s\n' % (len(ref)))

if __name__ == '__main__':
    main()
