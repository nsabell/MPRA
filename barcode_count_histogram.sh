#!/bin/bash

cd /home/nsabell/scratch/mpra

gunzip -c fastq_pairMerged/1KG_pairMerged_processed.fastq.gz | \
grep "@NS500418:" | cut -d : -f 8 | sort | uniq -c | \
sed 's/^ *//' | sed 's/ /\t/' > 1KG_pairMerged_processed_barcodeUniqFreq.txt &

gunzip -c fastq_pairMerged/GTEX_pairMerged_processed.fastq.gz | \
grep "@NS500418:" | cut -d : -f 8 | sort | uniq -c | \
sed 's/^ *//' | sed 's/ /\t/' > GTEX_pairMerged_processed_barcodeUniqFreq.txt &

wait