#!/usr/bin/env python

# Nathan Abell
# Montgomery Lab
# Stanford University
# Durga

# Recover oligo->barcode mappings

import sys

with open(sys.argv[1], 'r') as file, open(sys.argv[2], 'w') as outfile:
	lastBarcode = ""
	lastLine = ""
	skip = 0

	for line in file:

		data = line.strip().split()

		if lastBarcode == data[0] and skip == 1:
			continue
		
		if lastBarcode == data[0] and skip == 0:
			skip =1
			continue

		if lastBarcode != data[0] and skip == 1:
			lastBarcode = data[0]
			lastLine = line
			skip = 0
			continue

		if lastBarcode != data[0] and skip == 0:
			outfile.write(lastLine)
			lastBarcode = data[0]
			lastLine = line
			continue

	if lastBarcode != data[0] and skip == 0:
		outfile.write(lastLine)

