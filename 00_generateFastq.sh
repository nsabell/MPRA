#!/bin/bash

# Nathan Abell
# Montgomery Lab
# Stanford University
# SCG4

# ./scripts/00_generateFastq.sh ./170323_NS500418_0566_AHYNWJBGXY ./fastq #./fastqc

RAW=$1
FASTQ=$2
FASTQC=$3

module load bcl2fastq2
module load fastqc

cd $RAW
mkdir $FASTQ
bcl2fastq -o $FASTQ > bcl2fastq2.log 2>&1

cd $FASTQ
for i in {"GTEX","1KG"}; do
	R1="$i*R1";
	R2="$i*R2";
	
	cat $R1 > ${i}_R1.fastq.gz &
	cat $R2 > ${i}_R2.fastq.gz &
	wait;
	
	fastqc ${i}_R1.fastq.gz > ${i}_R1.fastqc.log 2>&1 &
	fastqc ${i}_R2.fastq.gz > ${i}_R2.fastqc.log 2>&1 &
	wait;

	mv *fastqc* $FASTQC;
done

