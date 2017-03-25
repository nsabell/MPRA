#!/bin/bash

# This is called by select_mpra_variants.py and not by anything else directly.
# This executes a join command on 

VCF=$1
VARIANTS=$2
OUTFILE=`basename $VARIANTS | cut -d . -f 1`

# Each variant in the original VCF and in the mora_test_variants.txt file is identified by a underscore separated chrom_pos, and then joined and formatted to produce the varInfo format
# Produces a file called mpra_test_variants.varInfo.txt
join <(awk '{print $1 "_" $2 "\t" $0}' $VARIANTS | sort -k1,1) <(cat $VCF | awk '{if ($1~/^[0-9]+$/) print $1 "_" $2 "\t" $3 "\t" $4 "\t" $5}' | sort -k1,1) | awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' > $OUTFILE.varInfo.txt


