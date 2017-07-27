# mpra

## MPRA Design from eQTL Summary Statistics and Population Genetic Data

## Random Barcode Data Analysis by Paired-End Sequencing

### Stage and Demultiplex

mkdir bin fastq fastqc logs raw
mkdir fastq/run1 fastq/run2 fastqc/run1 fastqc/run2

./scripts/00_generateFastq.sh ./170323_NS500418_0566_AHYNWJBGXY ./fastq ./fastqc
./scripts/00_generateFastq.sh ./170718_NS500418_0640_AHF7W5BGX2 ./fastq ./fastqc

mv ./170323_NS500418_0566_AHYNWJBGXY/fastq/*fastq.gz fastq/run1
mv ./170718_NS500418_0640_AHF7W5BGX2/fastq/*fastq/gz fastq/run2
mv ./170323_NS500418_0566_AHYNWJBGXY/fastqc/*fastqc* fastqc/run1
mv ./170718_NS500418_0640_AHF7W5BGX2/fastqc/*fastqc* fastqc/run2
mv ./170323_NS500418_0566_AHYNWJBGXY/fastq/DemultiplexReports
mv ./170718_NS500418_0640_AHF7W5BGX2/fastq/DemultiplexReports

tar cvfz 170323_NS500418_0566_AHYNWJBGXY.tar.gz 170323_NS500418_0566_AHYNWJBGXY  
tar cvfz 170718_NS500418_0640_AHF7W5BGX2.tar.gz 170718_NS500418_0640_AHF7W5BGX2

### Merge with FLASH and Get Unique Barcodes

01_readPairProcessing.sh <R1.fastq.gz> <R2.fastq.gz> <out_pfx>

### Pre-process For Alignment

### Create Oligo Reference

### Map Against Oligo Reference

Identification of Expression Modulating Variants

Nathan Abell

Montgomery Lab

Stanford University
