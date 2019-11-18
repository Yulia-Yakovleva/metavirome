from Bio.Blast import NCBIXML
from datetime import datetime

import csv


for i in range(0, 208):
    sample_list = []
    print(f'{datetime.now()}: {i} of 207')
    blast_records = NCBIXML.parse(open(f'out_blast_{i}.xml', 'r'))

    for blast_record in blast_records:
        query_id = blast_record.query.split(';')[0]
        for alignment in blast_record.alignments:
            acc = alignment.accession
            for hsp in alignment.hsps:
                cov = hsp.align_length / blast_record.query_length
                ident = hsp.identities / hsp.align_length
                evalue = hsp.expect
                sample_list.append([query_id, acc, cov, ident, evalue])

    with open(f'parse_out_{i}.csv', 'w') as f:
        f.write("Query_Name," "Hit_id," + "Query_cover," + "Identity," + "Evalue" + "\n")
        wr = csv.writer(f, dialect='excel')
        wr.writerows(sample_list)
