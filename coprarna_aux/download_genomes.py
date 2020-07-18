import os
import argparse
from Bio import Entrez
from urllib.request import HTTPError

information = """
The tool allows to download organism in gbk format and then zip it.
It downloads in current directory by default. You can also download at specified path. 

Usage:
Input:
* Organism id (Accession No.) 
* Organism destination folder (default: current directory) 
Output:
* Organism in zipped (.gz) format.

Example1:
- python download_genomes.py -g NC_000913 -p
"""


def download(organism, organism_path):

    fasta_file = organism
    Entrez.email = 'whatever@mail.com'

    f_name = organism_path + fasta_file + ".gb.gz"
    if(os.path.isfile(f_name)):
        return

    # accession id works, returns genome in fasta format, looks in the 'nucleotide' database:
    try:
        handle = Entrez.efetch(db='nucleotide', id=fasta_file,  rettype='gbwithparts', retmode="text", seq_start=1)  #,  seq_start=1 The genome with the accession number is fetched from ncbi server and saved in 'handle' variable
        # store locally:
        file_name = organism_path + fasta_file + ".gb"
        # print(file_name)
        local_file=open(file_name, 'w')
        local_file.write(handle.read()) # write the genome to a file
        handle.close()
        local_file.close()
        print(fasta_file + " succesfully downloaded")
    except HTTPError:
        print(fasta_file, "Error received in retrieving the fasta file with the given Accession Number.\nPlease check the Accession Number, or\nPlease check your internet connection or The 'NCBI' web server is down.\n")
        with open("organims_not_downloaded.txt", 'a') as wr:
            wr.write(fasta_file + '\n')



def zip_file(organism, organism_path):
    file_name = organism_path + organism.strip() + ".gb"
    if (os.stat(file_name).st_size != 0):
        cmd = 'gzip -9 ' + file_name
        # print(cmd)
        os.system(cmd)
    else:
        download(organism, organism_path)


def main():

    cur_dirr = os.popen("pwd").read()
    # print("Current Directory: ", cur_dirr)

    parser = argparse.ArgumentParser(prog="Script",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=information)
    parser.add_argument('-g', '--genome', type=str, nargs=1,
                        help='please enter the complete genome file name only')
    parser.add_argument('-p', '--gbk_path', type=str, nargs=1,
                        default=cur_dirr,
                        help='please enter the complete gbk file path')


    args = parser.parse_args()
    # parser.print_help()

    # path of the dna(s) file
    organism = args.genome[0]
    gbk_path = args.gbk_path[0]
    
    # print(organism)
    # print(gbk_path)

    # if complete genome protein extraction to be done
    if args.genome is not None:
        download(organism, gbk_path)
    else:
        print("please enter the organism to be downloaded")

    zip_file(organism, gbk_path)

if __name__ == "__main__":
    main()
