import argparse

from Bio import SeqIO

parser = argparse.ArgumentParser(
    description='Filter protein fasta file by keeping only proteins, that presented in other nucleotide fasta file.')

parser.add_argument('--nucleotide_fasta', '-n', required=True, help='Path to nucleotide fasta file')
parser.add_argument('--protein_fasta', '-p', required=True, help='Path to protein fasta file')
parser.add_argument('--output', '-o', required=True, help='Path to output protein fasta file')
args = parser.parse_args()
d = vars(args)

if __name__ == "__main__":
    nucleotide_contigs = SeqIO.parse(open(d['nucleotide_fasta']), 'fasta')
    protein_contigs = SeqIO.parse(open(d['protein_fasta']), 'fasta')
    filtered_protein_contigs = []

    allowed_contig_names = set([contig.name.split(';')[0] for contig in nucleotide_contigs])

    for protein_contig in protein_contigs:
        if protein_contig.name.split('_')[0] in allowed_contig_names:
            filtered_protein_contigs.append(protein_contig)

    SeqIO.write(filtered_protein_contigs, open(d['output'], mode='w'), 'fasta')
