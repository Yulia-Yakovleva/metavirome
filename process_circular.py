import sys
import glob
import os
import logging

from Bio import SeqIO

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

directory = 'seawater_metagenome/'


def is_circular(seq):
    for kval in range(50, 200):
        start = seq[:kval]
        end = seq[-kval:]
        if start == end:
            return True
    return False


if __name__ == "__main__":
    for i, subdir in enumerate(os.listdir(directory)):
        logging.info(f'{i} of {len(os.listdir(directory))}: {subdir}')

        subdir_path = directory + subdir
        if os.path.isfile(subdir_path):
            continue
        for file_name in os.listdir(subdir_path):
            file_path = subdir_path + '/' + file_name
            if file_name.endswith(".fna") and not file_name.endswith("_circular_only.fna") and not any(map(
                    lambda x: x == file_name.replace('.fna', '_circular_only.fna'), # check if it's file is already processed
                    os.listdir(subdir_path))):
                logging.info(f'<processing file>: {file_name}')

                contigs = SeqIO.parse(open(file_path, mode='r'), 'fasta')
                circular_contigs = []

                for contig in contigs:
                    seq = contig.seq
                    if len(seq) >= 500 and is_circular(seq):
                        print(file_name, contig.name)
                        circular_contigs.append(contig)

                SeqIO.write(circular_contigs, open(file_path.replace('.fna', '_circular_only.fna'), mode='w'), 'fasta')