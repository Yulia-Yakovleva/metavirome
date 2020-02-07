import os
import argparse

from Bio import SeqIO

parser = argparse.ArgumentParser(
    description='Parse table with `hmmsearch` results and extract proteins to output folder.')
parser.add_argument('--domtblout', '-d', required=True, help='Domtblout file')
parser.add_argument('--proteins', '-p', required=True, help='Proteins fasta')
parser.add_argument('--output', '-o', required=True, help='Output folfer name')
d = vars(parser.parse_args())
domtblout_file, proteins_file, outdir_name = d['domtblout'], d['proteins'], d['output']

# Get dict "HMM:proteins"
hmm_out = {}
with open(domtblout_file, "r") as domtblout:
    hmm_hits = domtblout.readlines()
hmm_hits = [i.split() for i in hmm_hits]

models_dict = {}

bitscore_thr = 20
eval_thr = 0.001

for i in hmm_hits:
    if i[0][0] != "#":
        if i[3] not in models_dict:
            models_dict[i[3]] = []
        if float(i[6]) < eval_thr and float(i[7]) > bitscore_thr:
            models_dict[i[3]].append(i[0])

# Get proteins
records = [s for s in SeqIO.parse(proteins_file, "fasta")]

# Create output directory
if not os.path.exists(outdir_name):
    os.mkdir(outdir_name)

# Write one fasta file per hmm
for i in models_dict:
    temp = []
    for rec in records:
        if rec.id in models_dict[i]:
            temp.append(rec)
    SeqIO.write(temp, outdir_name + "/" + i + ".fasta", "fasta")
