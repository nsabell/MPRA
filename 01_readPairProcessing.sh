#!/bin/bash

# Nathan Abell
# Montgomery Lab
# Stanford University
# SCG4

# Usage: 01_readPairProcessing.sh <R1.fastq.gz> <R2.fastq.gz> <out_pfx>

R1=$1
R2=$2
OUTPFX=$3
OUTDIR=`dirname $3`

module load fastx_toolkit
cd /home/nsabell/scratch/mpra
export FLASH=$PATH:/home/nsabell/scratch/mpra/bin/FLASH-1.2.11/flash
mkdir -p $OUTDIR

# Merge 1KG samples with Flash
$FLASH -z -r 150 -f 223 -s 10 fastq/run1/1KG_R1.fastq.gz fastq/run1/1KG_R2.fastq.gz
mv out.extendedFrags.fastq.gz fastq_pairMerged/run1/1KG_pairMerged.fastq.gz
mv out.notCombined_1.fastq.gz fastq_pairMerged/run1/1KG_R1_nonmerged.fastq.gz
mv out.notCombined_2.fastq.gz fastq_pairMerged/run1/1KG_R2_nonmerged.fastq.gz
mv out.hist fastq_pairMerged/run1/1KG_out.hist
mv out.histogram fastq_pairMerged/run1/1KG_out.histogram

mv *.o* *.e* logs

## Extract barcodes and collapse to get frequency
gunzip -c fastq_pairMerged/run1/1KG_pairMerged.fastq.gz | fastx_trimmer -l 20 | fastq_to_fasta -n -r | fastx_collapser | gzip > uniqBarcodes/1KG_uniqBarcodes.fastq.gz
gunzip -c fastq_pairMerged/run2/1KG_pairMerged.fastq.gz | fastx_trimmer -l 20 | fastq_to_fasta -n -r | fastx_collapser | gzip > uniqBarcodes/1KG_uniqBarcodes_combined.fastq.gz
gunzip -c fastq_pairMerged/GTEX_pairMerged.fastq.gz | fastx_trimmer -l 20 | fastq_to_fasta -n -r | fastx_collapser | gzip > uniqBarcodes/GTEX_uniqBarcodes.fastq.gz





