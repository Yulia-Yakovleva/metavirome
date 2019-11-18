from Bio import Entrez
from Bio import SeqIO
Entrez.email = "yakovleva.spbu@gmail.com"
import pandas as pd


my_csv = pd.read_csv('/home/yulia/metavirome/marin_and_seawater/blast_out/parse_out_110.csv')
column = my_csv.Hit_id
with open('/home/yulia/metavirome/marin_and_seawater/blast_out/parse_out_110.fasta', 'w') as out_file:
    for hit_id in set(column):
        # print(hit_id) # работаем с каждым айдишником хита
        with Entrez.efetch(db="protein", rettype="fasta", retmode="fasta", id=hit_id) as handle:
            seq_record = SeqIO.read(handle, "fasta")
        # print(seq_record.seq)
        SeqIO.write(seq_record, out_file, "fasta") 
