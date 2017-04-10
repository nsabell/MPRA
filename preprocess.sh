#!/bin/bash

module load python/2.7.6

./scripts/preprocess_merged_pairs.py ./fastq_pairMerged/1KG_pairMerged.fastq.gz &
./scripts/preprocess_merged_pairs.py ./fastq_pairMerged/GTEX_pairMerged.fastq.gz &

wait

