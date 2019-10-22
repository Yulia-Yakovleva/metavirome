#Usage: python domtblout_to_prot.py <domtblout_file> <proteins_fasta> <output_folder>


import sys
import os

domtblout_file = sys.argv[1]
proteins_file = sys.argv[2]
outdir_name = sys.argv[3]


# Get dict "HMM:proteins"
hmm_out={}
with open (domtblout_file, "r") as domtblout:
    hmm_hits = domtblout.readlines()
hmm_hits = [i.split() for i in hmm_hits]

models_dict = {}    
    
bitscore_thr = 20 
eval_thr = 0.001

for i in hmm_hits:
    if i[0][0]!="#":
        if i[3] not in models_dict:
            models_dict[i[3]] = []
        if float(i[6]) < eval_thr and float(i[7]) > bitscore_thr:
            models_dict[i[3]].append(i[0])
  


# Get proteins
records=[]
from Bio import SeqIO
for seq_record in SeqIO.parse(proteins_file, "fasta"):
    records.append(seq_record)


# Create output directory
if not os.path.exists (outdir_name):
  os.mkdir(outdir_name)

# Write one fasta file per hmm
for i in models_dict:
    temp = []
    for rec in records:        
            if rec.id in models_dict[i]:
                temp.append(rec)
    SeqIO.write(temp, outdir_name +"/" + i+".fasta", "fasta")


