# /usr/bin/python

from Bio import Entrez

import argparse

# argument parsing
parser = argparse.ArgumentParser(description='Parse information from NCBI Assembly database by given term and output it to the resulting text file.')
parser.add_argument('--term', '-t', required=True, help='Term, string value')
parser.add_argument('--output', '-o', required=True, help='Output file name')
d = vars(parser.parse_args())
term, output_file = d['term'], d['output']

# NCBI search
Entrez.email = "mike.rayko@gmail.com"
handle = Entrez.esearch(db="assembly", term=term, retmax=10500)
record = Entrez.read(handle)
print("got assembly ids")

id_list = record["IdList"]
print(len(id_list))
i = 0

spl = []
spl.append(id_list[:5000])
spl.append(id_list[5000:])
sample_list = []

# Collecting results from response
for id_l in spl:
    joined = ",".join(id_l)
    handle = Entrez.esummary(db="assembly", id=joined, retmax=10500)
    entry = Entrez.read(handle)
    elems = entry["DocumentSummarySet"]["DocumentSummary"]
    print("got summaries")
    print(len(elems))
    for elem in elems:
        i += 1
        if i % 500 == 0:
            print
            i
        taxid = elem["Taxid"]
        gbid = elem["GbUid"]
        acc = elem["AssemblyAccession"]
        taxid = elem["Taxid"]
        orgn = elem["Organism"]
        sptaxid = elem["SpeciesTaxid"]
        spname = elem["SpeciesName"]
        size = elem["Meta"].split("</Stat>")[14].split(">")[-1]

        sample_list.append([gbid, acc, taxid, orgn, sptaxid, spname, size])

# Writing results
with open(output_file, 'w') as f:
    f.write("# GbID,Accession,TaxID,Organism,SpeciesTaxID,SpeciesName,Size" + "\n")
    for item in sample_list:
        f.write("%s\n" % ",".join(item))
