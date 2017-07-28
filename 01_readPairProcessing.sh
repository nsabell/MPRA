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
export FLASH=/home/nsabell/scratch/mpra/bin/FLASH-1.2.11/flash
mkdir -p $OUTDIR

# Merge 1KG samples with Flash
$FLASH -o $OUTPFX -z -r 150 -f 223 -s 10 $R1 $R2
mv *.o* *.e* logs

## Extract barcodes and collapse to get frequency
gunzip -c ${OUTPFX}_pairMerged.fastq.gz | fastx_trimmer -l 20 | fastq_to_fasta -n -r | fastx_collapser | gzip > ${OUTPFX}_uniqBarcodes.fastq.gz


