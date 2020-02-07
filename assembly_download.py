import ftplib
import shutil
import os.path
import gzip
import logging
import argparse

from Bio import Entrez
from ftplib import FTP

# Parsing arguments
parser = argparse.ArgumentParser(description='Download assemblies from NCBI Assembly database by given term by accessing FTP folder and output it to the output folder.')
parser.add_argument('--term', '-t', required=True, help='Term, string value')
parser.add_argument('--output', '-o', required=True, help='Output folfer name')
d = vars(parser.parse_args())
term, output_dir = d['term'], d['output']

# Logger initiantion
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

# NCBI search
Entrez.email = 'a.zabelkin@itmo.ru'
general_output_dir = output_dir + '/'

handle = Entrez.esearch(db='assembly', term=term, retmax=10500)
record = Entrez.read(handle)
print('got assembly ids')

id_list = record['IdList']
print('len:', len(id_list))

# NCBI request for extra info
joined = ','.join(id_list)
handle = Entrez.esummary(db='assembly', id=joined, retmax=10500)
entry = Entrez.read(handle)
elems = entry['DocumentSummarySet']['DocumentSummary']

# FTP connection
ftp = FTP('ftp.ncbi.nlm.nih.gov')
ftp.login()

for i, elem in enumerate(elems):
    ftp_url = elem['FtpPath_GenBank']
    if ftp_url == '': continue

    dir_path = ftp_url.replace('ftp://ftp.ncbi.nlm.nih.gov', '')
    logging.info(f'{i} of {len(elems)}: {dir_path}')
    ftp.cwd(dir_path)

    # For every file in folder
    for file in ftp.nlst():
        output_dir = general_output_dir + elem['AssemblyAccession'] + '/'
        output_file = output_dir + file

        # Creating directory
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # If it's not truly nucleotide fasta
        if '.fna' in file \
                and not 'rna_from' in file \
                and not 'cds_from' in file \
                and not os.path.isfile(output_file):
            # Downloading
            logging.info(f'<downloading file>: {file}')
            try:
                with open(output_file, 'wb') as local_file:
                    ftp.retrbinary('RETR ' + file, local_file.write)
            except ftplib.error_perm:
                pass

ftp.quit()
