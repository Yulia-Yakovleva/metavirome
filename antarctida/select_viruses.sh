#!/bin/sh
hmmsearch --domtblout $1.out --incT 30 --cpu 8 ../viral_hmms/selected_profiles_Sept2019.hmm "/run/media/my/New Volume/temp/antarctica_metagenomes/"$1.fasta
python ../domtblout_to_prot.py $1.out "/run/media/my/New Volume/temp/antarctica_metagenomes/"$1.fasta $1

