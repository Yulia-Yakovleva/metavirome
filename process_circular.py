import os
import logging
import gzip
import shutil
import argparse

from Bio import SeqIO

# logger initiation
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

# argument parsing
parser = argparse.ArgumentParser(description='Filter circular contigs in folder produced by `assembly_download.py` script.')
parser.add_argument('--folder', '-d', required=True, help='Folder with assemblies')
parser.add_argument('--remove', '-o', required=True, type=bool, help='Remove archive files')
d = vars(parser.parse_args())
directory, remove_archive_files = d['folder'], d['remove']

# Classic checking for circular
def is_circular(seq):
    for kval in range(50, 200):
        start = seq[:kval]
        end = seq[-kval:]
        if start == end:
            return True
    return False


# Same, but Using Knuth–Morris–Pratt algorithm
def is_circular_knp(seq, min_match=50, max_match=200):
    def prefix(s):
        v = [0] * len(s)
        for i in range(1, len(s)):
            k = v[i - 1]
            while k > 0 and s[k] != s[i]:
                k = v[k - 1]
            if s[k] == s[i]:
                k += 1
            v[i] = k
        return v

    s = seq[:max_match] + "$" + seq[-max_match:]
    ps = prefix(s)

    return ps[-1] >= min_match


if __name__ == "__main__":
    # For every subdirectory
    for i, subdir in enumerate(os.listdir(directory)):
        logging.info(f'{i} of {len(os.listdir(directory))}: {subdir}')

        subdir_path = directory + subdir
        if os.path.isfile(subdir_path): continue
        # For every file in directory
        for file_name in os.listdir(subdir_path):
            file_path = subdir_path + '/' + file_name
            circular_file_name = file_name.replace('.fna', '_circular_only.fna').replace('.gz', '')
            circular_file_path = subdir_path + '/' + circular_file_name

            if file_name.endswith(".fna.gz") and not any(map(
                    lambda x: x == circular_file_name, # check if it's file is already processed
                    os.listdir(subdir_path))):

                # Unarchiving file
                logging.info(f'<unarchiving file>: {file_name}')
                with gzip.open(file_path, 'rb') as f_in:
                    with open(file_path[:-3], 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                contigs = SeqIO.parse(open(file_path[:-3]), 'fasta')

                # Processing file
                logging.info(f'<processing file>: {file_name}')
                circular_contigs = []
                for contig in contigs:
                    seq = contig.seq
                    if len(seq) >= 500 and is_circular_knp(seq):
                        print(file_name, contig.name)
                        circular_contigs.append(contig)

                SeqIO.write(circular_contigs, open(circular_file_path, mode='w'), 'fasta')

                # Remove unarchived file
                os.remove(file_path[:-3])
                logging.info(f'<unarchived file removed>: {file_path[:-3]}')
                if remove_archive_files:
                    os.remove(file_path)
                    logging.info(f'<archive file removed>: {file_path}')