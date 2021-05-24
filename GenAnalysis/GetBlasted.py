### Script which extracts sequences according to blast search
### example to run: python3 GetBlasted.py <reference MG1655 sequences - fasta file> <path to blast results> <path to env. strains' genomes in fasta> <path to output directory> <number of which alignment length can be smaller than MG1655 reference to be included>

import sys
import glob
from Bio.Seq import reverse_complement
from Bio import SeqIO, SeqFeature
from tqdm import tqdm

### function to determine number of lines in a file
def get_num_lines(file_path):
    query = open(file_path, 'r')
    lines = 0
    for line in query:
        lines += 1
    return lines

### main function
def main():

    ### save arguments
    k12_proms = sys.argv[1]
    blast_dir = sys.argv[2]
    strains_dir = sys.argv[3]
    out_dir = sys.argv[4]
    diff = sys.argv[5]

    ### read all files with blast results
    blast_results = glob.glob(blast_dir + '*.txt')

    print("\nExtracting blasted sequences.")
    pbar = tqdm(blast_results)
    ### go throught all strain files with blast results
    for strain in pbar:
        str = strain.split("/")[len(strain.split("/"))-1]
        str_only = str.split("_")[0]
        print(str_only)
        pbar.set_description("Processing %s" % str)
        with open(strain, 'r') as file:
            ### make a separate output file for each strain
            prom_out = open(out_dir + '{}_-{}.fasta'.format(str_only, diff), 'w')
            ### loop through all lines in strain's blast results
            for row in file:
                ### save important variables from each blast result
                prom_idf = row.split("\t")[0]
                node_idf = row.split("\t")[1]
                align_length = row.split("\t")[3]
                sstart = int(row.split("\t")[8])
                send = int(row.split("\t")[9])
                evalue = row.split("\t")[10]
                ### now loop through all extracted MG1655 reference sequences
                with open(k12_proms, 'r') as proms:
                    for line in proms:
                        ### save promoter name if line starts with ">"
                        if line.startswith(">"):
                            prom_idp = line.split(">")[1].split("\n")[0]
                        ### if not
                        else:
                            ### save length of the reference sequence
                            prom_length = len(line)
                            ### continue only if promoter names for reference and env. strain are identical
                            ### and at the same time the length of the blasted sequence is not shorter then 'diff'bp as compared to the reference
                            if all([prom_idf == prom_idp, int(align_length) >= prom_length-int(diff)]):
                                ### loop through all sequences of that particular strain in its fasta file
                                for fsa in SeqIO.parse(strains_dir + str_only + "/" + str_only + '.fsa', "fasta"):
                                    ### find the sequence where blast found a hit with the MG1655 reference
                                    if node_idf in repr(fsa.id):
                                        ### and extract the sequence
                                        prom_out.write(">" + fsa.id + "|" + prom_idp + "|alignment_length:_" + align_length + "|evalue:_" + evalue + "|strand:_")
                                        if sstart < send:
                                            prom_out.write("+\n%s\n" % fsa.seq[sstart-1:send])
                                        else:
                                            prom_out.write("-\n%s\n" % reverse_complement(fsa.seq[send-1:sstart]))

    ### read all files from previous step
    blast_extracted = glob.glob(out_dir + '*{}.fasta'.format(diff))

    print("\nCreating separated file for each promoter.")
    pbar = tqdm(SeqIO.parse(k12_proms, "fasta"), total=get_num_lines(k12_proms)//2)
    ### loop through all lines in file with extracted MG1655 reference sequences
    for entry in pbar:
        ### save promoter name only
        promoter = entry.id.split("_")[3]
        pbar.set_description("Processing %s" % promoter)
        ### create a separate file for each promoter name
        out_by_prom = open(out_dir + 'by_promoter/' + '{}_-{}.fasta'.format(promoter, diff), 'w')
        ### write MG1655 annotation and sequence into it
        out_by_prom.write(">%s\n%s\n" % (entry.id, entry.seq))
        ### loop through all sequence extracted in previous loop for all strains
        for file in blast_extracted:
            file_only = file.split("/")[len(file.split("/"))-1]
            str_id = file_only.split("_")[0]
            for line in SeqIO.parse(file, "fasta"):
                ### if there is such homologous sequence present in that file, add it to the output file
                if "_{}_".format(promoter) in repr(line.id):
                    out_by_prom.write(">%s|%s\n%s\n" % (str_id, line.id.split("\n")[0], line.seq))


if __name__ == '__main__':
    main()
