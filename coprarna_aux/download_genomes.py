import argparse
import os
import wget
import sys
import signal
import urllib.request
import threading
from time import sleep
import shutil

# absolute path of this script
scriptPath = os.path.dirname(os.path.realpath(__file__))

# GLOBAL VARIABLES
genbank_summary_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/" +\
                    "prokaryotes.txt"
ncbi_folder_name = scriptPath + "/myNCBI-LOOKUP"
ncbi_master_table = "NCBI-MASTER.table"

# disabling ctrl z option
signal.signal(signal.SIGTSTP, signal.SIG_IGN)


def get_data_from_ncbi(ftp_lnk, user_acc, flag, download_path):
    try:
        if flag:
            organism_path = str(download_path) + user_acc.strip().split(".")[0] + ".gbk.gz"
            os.system("wget -c --retry-connrefused -nv --show-progress --continue\
                      --read-timeout=20 --tries=40 --wait=10 --timeout=15 " +
                      str(ftp_lnk.strip()) + " -O " + str(organism_path))
            print()
            os.system("gunzip " + organism_path)
            with open(organism_path[:-3], 'r') as f:
                file_open_flag = False
                output_file = ""
                for line in f:
                    if "LOCUS" in line:
                        file_open_flag = True
                        user_acc = line.strip().split()[1]
                        output_file = str(download_path) + user_acc.strip().split(".")[0] + ".gb"
                        wr = open(output_file, 'a+')
                        wr.write(line)
                    elif not line.startswith("//") and file_open_flag:
                        wr.write(line)
                    elif line.startswith("//") and file_open_flag:
                        file_open_flag = False
                        wr.write(line)
                        wr.close()
                        zip_file(output_file)

            os.system("rm -f " + organism_path[:-3])

        else:
            org_link = '"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db'\
                       '=nucleotide&id=' + user_acc + '&rettype=gbwithparts&retmode=text"'
            # output file name with location
            output_file = str(download_path) + "/" + user_acc.strip().split(".")[0] + ".gb"
            os.system("wget -c --retry-connrefused -nv --show-progress --continue\
                      --read-timeout=20 --tries=40 --wait=10 --timeout=15 " +
                      str(org_link) + " -O " + str(output_file))
            zip_file(output_file)

    except Exception:
        sys.stderr.write("***FTP/HTTP ERROR: " + user_acc + " CANNOT BE DOWNLOADED ***")
        shutil.rmtree(output_file)


def get_data(user_acc, organism_col, ftp_lnk, download_path):
    try:
        urllib.request.urlopen(ftp_lnk)
        if user_acc in organism_col:
            flag = True
            get_data_from_ncbi(ftp_lnk, user_acc, flag, download_path)

    except Exception:
            flag = False
            # chr_ref =  user_acc
            organism_col = user_acc
            ftp_lnk = 'NA'
            # chr_flag = "NA"
            get_data_from_ncbi(ftp_lnk,  user_acc, flag, download_path)


def zip_file(organism):
    file_name = organism.strip()
    outpath = file_name + ".gz"
    if not os.path.exists(outpath):
        if (os.stat(file_name).st_size != 0):
            cmd = 'gzip -9 ' + file_name
            os.system(cmd)
        else:
            find_link(organism)


def process_NCBI_lookup(in_file, directory, out_file):
    """ Downloads prokayrotes.txt from a ftp link and then
        extracts the accession numbers and ftp link to create
        a new file for organisms' lookup"""
    result = list()
    handle = open(in_file, "r", encoding='ascii', errors="replace")
    for line in handle:
        line = line.rstrip()
        tmp_arr = line.split("\t")
        # choosing only the entries with chromosomes.
        if "chromosome" in line and not line.startswith("#") and (len(tmp_arr) == 23):
            try:
                tmp_identifiers_arr = tmp_arr[8].split("; ")
                file_identifiers = ""
                chromosome = []
                plasmid = []
                for line in tmp_identifiers_arr:
                    # splitting the replicons to later write to a file
                    replicons = line.strip().split(":")
                    # if list item has chromosome
                    if "chromosome" in replicons[0].strip().split()[0]:
                        chromosome.append(replicons[1].strip())
                    # if list item has plasmid
                    if "plasmid" in replicons[0].strip().split()[0]:
                        plasmid.append(replicons[1].strip())
                # if no plasmid or chromosome entry, then put "-"
                if len(plasmid) == 0:
                    plasmid.append("-")
                if len(chromosome) == 0:
                    chromosome.append("-")
                # prefix of file to download
                missing_part = tmp_arr[20].split("/")[-1]
                file_identifiers = ','.join(chromosome) + "\t" + ','.join(plasmid) +\
                                   "\t" + str(tmp_arr[20]) + "/" + str(missing_part) +\
                                   "_genomic.gbff.gz" + "\n"
                result.append(file_identifiers)
            except Exception:
                print("Prokayrotes file is missing at FTP link !!!")

    # write output
    out_path = str(directory) + "/" + str(out_file)
    handle = open(out_path, "w")
    for entry in result:
        handle.write(entry)
    handle.close()
    # clean prokayrotes.txt file
    os.remove(in_file)


def find_link(input_param):
    # User starts a search - check if the ACC is in the database or not.
    # If not, download the genbank from ncbi
    input_param = input_param.strip().split(".")[0]
    handle = open(ncbi_file_path, "r")
    organism_col = ''
    ftp_lnk = ''
    for line in handle:
        line = line.rstrip()
        line_arr = line.split("\t")
        for i, val in enumerate(line_arr[:2]):
            if input_param in val:
                organism_col = val
                ftp_lnk = line_arr[2].strip()
                if input_param == val.strip().split("/")[0].strip().split(".")[0]:
                    ftp_lnk = ftp_lnk.replace('GCA', 'GCF')
    handle.close()

    if args.genomePath:
        download_path = str(args.genomePath) + "/"
    else:
        download_path = os.getcwd() + "/"
    if not os.path.isdir(download_path):
        os.system("mkdir " + str(download_path))
    file_path = str(download_path) + input_param.strip().split(".")[0] + ".gb.gz"
    if not os.path.exists(file_path):
        get_data(input_param, organism_col, ftp_lnk, download_path)

    # for threading to introduce delay
    sGlobal.release()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--accession", nargs="+", help="accession id\
                        e.g. CP021219.1", type=str, default="")
    parser.add_argument("-p", "--genomePath", help="Genome Path\
                        directory ", type=str, default="")
    parser.add_argument("-c", "--cores", help="no. of cores\
                        e.g 10", type=int, default="3")

    args = parser.parse_args()
    # ncbi_file_path = ''
    ncbi_file_path = str(ncbi_folder_name) + "/" + str(ncbi_master_table)

    # check if ncbi lookup table exist,
    # otherwise download the table from ncbi's ftp server
    if os.path.isdir(ncbi_folder_name):
        ncbi_master_table_path = str(ncbi_folder_name) + "/" + str(ncbi_master_table)
        if os.path.isfile(ncbi_master_table_path):
            pass
        else:
            # get data and reprocess the file
            prok_file = str(ncbi_folder_name) + "/prokaryotes.txt"
            file_path = ''
            if not os.path.isfile(prok_file):
                file_path = wget.download(
                            genbank_summary_url, out=ncbi_folder_name
                            )
            else:
                file_path = prok_file
            process_NCBI_lookup(file_path, ncbi_folder_name, ncbi_master_table)
            print()
    else:
        os.mkdir(ncbi_folder_name)
        # get data and reprocess the file
        file_path = wget.download(genbank_summary_url, out=ncbi_folder_name)
        process_NCBI_lookup(file_path, ncbi_folder_name, ncbi_master_table)
        print()

    global sGlobal
    threads = []
    sGlobal = threading.Semaphore(args.cores)
    if args.accession:
        for org in args.accession:
            sGlobal.acquire()
            t = threading.Thread(target=find_link, args=(org,))
            threads.append(t)
            t.start()
            sleep(0.34)
        for t in threads:
            t.join()
