#!/bin/bash

# Nathan Abell
# Montgomery Lab
# Stanford University
# SCG4

module load fastx_toolkit
module load bwa/0.7.15

cd /home/nsabell/scratch/mpra

cat reference/fasta/*1KG* reference/fasta/*sabeti* | sed '/--/d' | fastx_trimmer -f 16 | fastx_trimmer -t 15 > reference/fasta/1KG_all.fa
cat reference/fasta/*rare* | sed '/--/d' | fastx_trimmer -f 16 | fastx_trimmer -t 15 > reference/fasta/GTEX_all.fa

./scripts/collapse_reference_fasta.py reference/fasta/1KG_all.fa
.scripts/collapse_reference_fasta.py reference/fasta/GTEX_all.fa

mkdir -p reference/bwa/1KG reference/bwa/GTEX

bwa index -p reference/bwa/1KG/1KG reference/fasta/1KG_all_collapsed.fa
bwa index -p reference/bwa/GTEX/GTEX reference/fasta/GTEX_all_collapsed.fa


