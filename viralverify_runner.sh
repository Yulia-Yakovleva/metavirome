#!/usr/bin/env bash
folder='seawater_metagenome'

for file in `find $folder -name "*_genomic_circular_only.fna"`; do
    dir_path=`dirname "$file"`
    python ~/viralVerify/viralverify.py -f ./${file} -o ./${dir_path}/viralverify_output --hmm ~/viralVerify/nbc_hmms.hmm -t 4 -p
done