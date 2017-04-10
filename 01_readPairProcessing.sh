#!/bin/bash

# Nathan Abell
# Montgomery Lab
# Stanford University
# SCG4

export FLASH=$PATH:/home/nsabell/scratch/mpra/bin/FLASH-1.2.11/flash
cd /home/nsabell/scratch/mpra
mkdir fastq_pairMerged logs

# Merge 1KG samples with Flash
$FLASH -z -r 150 -f 223 -s 10 fastq/1KG_R1.fastq.gz fastq/1KG_R1.fastq.gz
mv out.extendedFrags.fastq.gz fastq_pairMerged/1KG_pairMerged.fastq.gz
mv out.notCombined_1.fastq.gz fastq_pairMerged/1KG_R1_nonmerged.fastq.gz
mv out.notCombined_2.fastq.gz fastq_pairMerged/1KG_R2_nonmerged.fastq.gz
mv out.hist fastq_pairMerged/1KG_out.hist
mv out.histogram fastq_pairMerged/1KG_out.histogram

# Merge GTEX samples with Flash
$FLASH -z -r 150 -f 223 -s 10 fastq/GTEX_R1.fastq.gz fastq/GTEX_R2.fastq.gz
mv out.extendedFrags.fastq.gz fastq_pairMerged/GTEX_pairMerged.fastq.gz
mv out.notCombined_1.fastq.gz fastq_pairMerged/GTEX_R1_nonmerged.fastq.gz
mv out.notCombined_2.fastq.gz fastq_pairMerged/GTEX_R2_nonmerged.fastq.gz
mv out.hist fastq_pairMerged/GTEX_out.hist
mv out.histogram fastq_pairMerged/GTEX_out.histogram

mv *.o* *.e* logs

## Extract barcodes and collapse to get frequency
gunzip -c fastq_pairMerged/1KG_pairMerged.fastq.gz | fastx_trimmer -l 20 | fastq_to_fasta -n -r | fastx_collapser | gzip > 1KG_uniqBarcodes.fastq.gz
gunzip -c fastq_pairMerged/GTEX_pairMerged.fastq.gz | fastx_trimmer -l 20 | fastq_to_fasta -n -r | fastx_collapser | gzip > GTEX_uniqBarcodes.fastq.gz





