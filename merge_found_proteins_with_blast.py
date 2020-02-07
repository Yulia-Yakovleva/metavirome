from Bio import SeqIO, Entrez
from datetime import datetime

import pandas as pd

Entrez.email = "yakovleva.spbu@gmail.com"

# This script just merge hmmsearch proteins with blast proteins and does not anything complex

sub_folder = 'marine_and_seawater_metagenome/jya_from_server/'
protein_names_file_path = sub_folder + 'protein_names.txt'


with open(protein_names_file_path) as f:
    protein_names = f.read().splitlines()

for ind, protein_name in enumerate(protein_names):
    print(f'{datetime.now()}: {ind}')

    blast_contigs_df = pd.read_csv(sub_folder + f'blast_out/parse_out_{ind}.csv')
    found_contigs = list(SeqIO.parse(open(sub_folder + 'domtblout_to_prot/' + protein_name), 'fasta'))

    joined = ','.join(blast_contigs_df['Hit_id'].unique())
    with Entrez.efetch(db="protein", rettype="fasta", retmode="fasta", id=joined) as handle:
        blast_contigs = SeqIO.parse(handle, "fasta")
        found_contigs += list(blast_contigs)

    SeqIO.write(found_contigs, open(
        sub_folder + 'found_proteins_with_blast_proteins/' + protein_name.replace('.fasta', '_with_blast_alignments.fasta'),
        mode='w'), 'fasta')
