#/usr/bin/python    

from Bio import Entrez
Entrez.email = "mike.rayko@gmail.com" 
handle = Entrez.esearch(db="assembly", term="metagenome",retmax = 10500) 
record = Entrez.read(handle)
print ("got assembly ids")

id_list = record["IdList"]
print (len(id_list))
i = 0

spl = []
spl.append(id_list[:5000])
spl.append(id_list [5000:])
sample_list=[]

for id_l in spl:
    joined=",".join(id_l)
    handle = Entrez.esummary(db="assembly", id=joined, retmax=10500)
    entry = Entrez.read(handle)
    elems = entry["DocumentSummarySet"]["DocumentSummary"]
    print ("got summaries")
    print (len(elems))
    for elem in elems:
        i += 1
        if i % 500 == 0:
            print(i)
        id = elem.attributes["uid"]
        taxid = elem["Taxid"]
        gbid = elem["GbUid"]
        acc = elem["AssemblyAccession"]
        taxid = elem["Taxid"]
        orgn = elem["Organism"]
        sptaxid = elem["SpeciesTaxid"]
        spname = elem["SpeciesName"]
        size = elem["Meta"].split("</Stat>")[14].split(">")[-1]

        sample_list.append([id, gbid, acc, taxid, orgn, sptaxid, spname, size])
    
with open('assembly_stats.txt', 'w') as f:
    f.write("# id, GbID,Accession,TaxID,Organism,SpeciesTaxID,SpeciesName,Size" + "\n")
    for item in sample_list:
        f.write("%s\n" % ",".join(item))
