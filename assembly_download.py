import ftplib
import shutil
import os.path
import gzip
import logging

from Bio import Entrez
from ftplib import FTP

term = 'soil metagenome'

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

Entrez.email = 'a.zabelkin@itmo.ru'
general_output_dir = term.replace(' ', '_') + '/'

handle = Entrez.esearch(db='assembly', term=term, retmax=10500)
record = Entrez.read(handle)
print('got assembly ids')

id_list = record['IdList']
print('len:', len(id_list))

joined = ','.join(id_list)
handle = Entrez.esummary(db='assembly', id=joined, retmax=10500)
entry = Entrez.read(handle)
elems = entry['DocumentSummarySet']['DocumentSummary']

ftp = FTP('ftp.ncbi.nlm.nih.gov')
ftp.login()

for i, elem in enumerate(elems):
    ftp_url = elem['FtpPath_GenBank']
    if ftp_url == '':
        continue

    dir_path = ftp_url.replace('ftp://ftp.ncbi.nlm.nih.gov', '')
    logging.info(f'{i} of {len(elems)}: {dir_path}')
    ftp.cwd(dir_path)

    files = ftp.nlst()
    for file in files:
        output_dir = general_output_dir + elem['AssemblyAccession'] + '/'
        output_file = output_dir + file

        # Creating directory
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        if '.fna' in file \
                and not 'rna_from' in file \
                and not 'cds_from' in file \
                and not 'spades' in file \
                and not os.path.isfile(output_file):
            # Downloading
            logging.info(f'<downloading file>: {file}')
            try:
                with open(output_file, 'wb') as local_file:
                    ftp.retrbinary('RETR ' + file, local_file.write)
            except ftplib.error_perm:
                pass

ftp.quit()
