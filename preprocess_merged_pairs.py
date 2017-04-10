#!/bin/env python

# Nathan Abell
# Montgomery Lab
# Stanford University
# SCG4

# preprocess_merged_pairs.py <input_fastq_from_flash>

import sys
import gzip

outname = (".").join(sys.argv[1].split('.')[:-2]) + "_processed.fastq.gz"

with gzip.open(outname, 'wb') as outfile, gzip.open(sys.argv[1], 'rb') as file:

	while True:

		header = file.readline()
		if not header:
			break
		header = header.strip().split(' ')[0]
		sequence = file.readline().strip()
		filter = file.readline()
		phred = file.readline().strip()

		outfile.write(header + ":" + sequence[0:20] + '\n' + sequence[58:-16] + '\n' + filter + phred[58:-16] + '\n')

