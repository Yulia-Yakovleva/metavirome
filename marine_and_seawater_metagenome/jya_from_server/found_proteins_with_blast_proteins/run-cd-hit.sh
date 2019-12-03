#!/usr/bin/env bash
for i in `find . -name "*.fasta" -type f`; do
    cd-hit -i $i -o "${i%.*}" -c 0.5 -n 3 -g 1 -aS 0.8 -T 4 -M 16000
done

rm *alignments
