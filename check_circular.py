import sys 
import os 
import string 
import re 
import subprocess 
import datetime 
import fastaparser 
from os.path import join 
from genericpath import isdir, exists


if __name__ == "__main__":
    file = sys.argv[1]
    contigs = fastaparser.read_fasta(file)
    count = []

    for contig in contigs:
        arr = contig[0].strip(';').split('_')
        for kval in range (200,50, -1):
            if kval >= len(contig[1]) or len(contig[1]) < 500:
                continue
            start = contig[1][:kval]
            end = contig[1][-kval:]
            if start == end:
                print (contig[0])
                print (contig[1])
                break   


