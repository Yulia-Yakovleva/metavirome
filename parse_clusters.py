import argparse

from glob import glob
from Bio import SeqIO

# Parsing arguments
parser = argparse.ArgumentParser(description='Parse clusters and filter fasta files by leaving only cluster centorid in given folder with protein .fasta files and them .clstr clusters file. '
                                             'New fasta centroids fasta file names have name equal name as original files, but with `_centroids` in the end.')
parser.add_argument('--folder', '-d', required=True, help='Working folder')
d = vars(parser.parse_args())
folder = d['folder']

# Parse cluster file
def parse_clusters_file(file_name):
    clusters = []

    with open(file_name, 'r') as f:
        lines = f.read().splitlines()

    for line in lines:
        if line[0] == ">":
            assert len(clusters) == int(line.split(' ')[1])
            clusters.append([])
        else:
            splitted = line.split(' ')
            name = splitted[1].split('...')[0].replace('>', '')
            if len(splitted) == 3:
                clusters[-1].insert(0, name)
            else:
                clusters[-1].append(name)

    return clusters

# Get centroids from fasta by given parsed clusters
def get_centroids_from_fasta(file_name, clusters):
    fasta_contigs = list(SeqIO.parse(open(file_name), 'fasta'))
    centroid_contigs = []

    for i, cluster in enumerate(clusters):

        is_b, is_m = False, False
        for clstr_contig in cluster:
            is_m |= ';size' in clstr_contig
            is_b |= not ';size' in clstr_contig

        centroid = cluster[0]

        # Here we can use binary search if we want better performance
        for fasta_contig in fasta_contigs:
            if centroid in fasta_contig.name:
                fasta_contig.id = f'{"B" if is_b else ""}{"M" if is_m else ""}_Centroid_{i}_Size_{len(cluster)}'
                fasta_contig.description = ''

                centroid_contigs.append(fasta_contig)
                break

    return centroid_contigs


clstr_files = glob(folder + '*.clstr')

for clstr_file in clstr_files:
    fasta_filename = clstr_file.replace('.clstr', '.fasta')
    fasta_centroids_filename = clstr_file.replace('.clstr', '_centroids.fasta')

    clusters = parse_clusters_file(clstr_file)
    centroid_contigs = get_centroids_from_fasta(fasta_filename, clusters)

    SeqIO.write(centroid_contigs, open(fasta_centroids_filename, mode='w'), 'fasta')
