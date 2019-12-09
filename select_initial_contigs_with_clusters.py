from parse_clusters import parse_clusters_file

from glob import glob
from pprint import pprint
from Bio import SeqIO
from collections import defaultdict


def collect_m_centroids(clusters, model, m_centroids):
    for i, cluster in enumerate(clusters):

        is_b, is_m = False, False
        for clstr_contig in cluster:
            is_m |= ';size' in clstr_contig
            is_b |= not ';size' in clstr_contig

        centroid = cluster[0].replace(';size', '')
        if is_m and not is_b:
            m_centroids[centroid].append(model + f'_Centroid_{i}_Size_{len(cluster)}')

all_contigs_filename = 'marine_and_seawater_metagenome/circular_merged.fna'
resulting_file = 'marine_and_seawater_metagenome/viral_annotated_circular_merged.fna'

folder = 'marine_and_seawater_metagenome/jya_from_server/found_proteins_with_blast_proteins/'
clstr_files = glob(folder + '*.clstr')
m_centroids = defaultdict(lambda: [])

for clstr_file in clstr_files:
    fasta_filename = clstr_file.replace('.clstr', '.fasta')
    fasta_centroids_filename = clstr_file.replace('.clstr', '_centroids.fasta')
    model = clstr_file.split('/')[-1].split('_')[0]

    clusters = parse_clusters_file(clstr_file)
    collect_m_centroids(clusters, model, m_centroids)

pprint(m_centroids)
all_circular_contigs = SeqIO.parse(open(all_contigs_filename), 'fasta')

m_contigs = []
for contig in all_circular_contigs:
    if contig.name in m_centroids.keys():
        id, desc = contig.description.split(' ', 1)

        models = ','.join(m_centroids[contig.name])
        contig.description = f'{id} Models:{models} {desc}'

        m_contigs.append(contig)

SeqIO.write(m_contigs, open(resulting_file, mode='w'), 'fasta')