# Black cat in a dark room: search for new viruses in metagenomes

###  Project goals
Viruses are the most abundant biological entities on Earth, but they are challenging in detecting, isolating, and classifying. 
The databases contain many assembled metagenomes, however, the diversity of viruses in assemblies is analyzed quite superficially.
The main goal is analyzing assembled metagenomes from NCBI Assembly database and preprocessed metagenomic raw reads with predicted proteins from MG-RAST database.

### Methods
We extracted сyclic sequences using a **Knuth-Morris-Pratt algorithm**. 
We determined whether the sequence belongs to the virus with **viralVerify** contig verification tool. 
From the contigs which were classified as viral we extracted structural genes like capsid genes and terminases using **hmmsearch**. 
The resulting sequences were aligned against nr **BLAST** database.
We clustered the resulting protein sequences with **CD-HIT**. 
The resulting centroid sequences were aligned using **MAFFT**. 
Then maximum likelihood trees were constructed by **RAxML**. 
Viruses contigs were annotated with **Prokka**. 
Annotations were visualised in **Geneious**.

## Requirements
All project-specific code was written with python language using biopython, argparse and pandas libraries.
You can use any platform with python and required libraries installed.

## Scripts description
Main scripts of project can be run with python and console arguments.

### assembly_download.py
This script downloads assemblies from NCBI Assembly database by given term by accessing FTP folder.

```bash
$ python assembly_download.py -h
usage: assembly_download.py [-h] --term TERM --output OUTPUT

Download assemblies from NCBI Assembly database by given term by accessing FTP
folder and output it to the output folder.

optional arguments:
  -h, --help            show this help message and exit
  --term TERM, -t TERM  Term, string value
  --output OUTPUT, -o OUTPUT
                        Output folfer name
```

### assembly_stats.py
This script parses information from NCBI Assembly database by given term and output it to the resulting text file.
```bash
$ python assembly_stats.py -h
usage: assembly_stats.py [-h] --term TERM --output OUTPUT

Parse information from NCBI Assembly database by given term and output it to
the resulting text file.

optional arguments:
  -h, --help            show this help message and exit
  --term TERM, -t TERM  Term, string value
  --output OUTPUT, -o OUTPUT
                        Output file name
```

### process_circular.py
This script filters circular contigs in folder produced by `assembly_download.py` script.

```bash
$ python process_circular.py -h
usage: process_circular.py [-h] --folder FOLDER --remove REMOVE

Filter circular contigs in folder produced by `assembly_download.py` script.

optional arguments:
  -h, --help            show this help message and exit
  --folder FOLDER, -d FOLDER
                        Folder with assemblies
  --remove REMOVE, -o REMOVE
                        Remove archive files
```

### domtblout_to_prot.py
This script parses table with `hmmsearch` results and extract proteins to output folder
```bash
$ python domtblout_to_prot.py -h
usage: domtblout_to_prot.py [-h] --domtblout DOMTBLOUT --proteins PROTEINS
                            --output OUTPUT

Parse table with `hmmsearch` results and extract proteins to output folder.

optional arguments:
  -h, --help            show this help message and exit
  --domtblout DOMTBLOUT, -d DOMTBLOUT
                        Domtblout file
  --proteins PROTEINS, -p PROTEINS
                        Proteins fasta
  --output OUTPUT, -o OUTPUT
                        Output folfer name
```

### parse_clusters.py
This script
```bash
$ python parse_clusters.py -h
usage: parse_clusters.py [-h] --folder FOLDER

Parse clusters and filter fasta files by leaving only cluster centorid in
given folder with protein .fasta files and them .clstr clusters file. New
fasta centroids fasta file names have name equal name as original files, but
with `_centroids` in the end.

optional arguments:
  -h, --help            show this help message and exit
  --folder FOLDER, -d FOLDER
                        Working folder
```