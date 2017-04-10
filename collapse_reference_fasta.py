#!/bin/env python

# Nathan Abell
# Montgomery Lab
# Stanford University
# SCG4

# collapse_reference_fasta.py <input_fasta>

import sys
from Bio import SeqIO

fasta_sequences = SeqIO.parse(open(sys.argv[1]), "fasta")
seqDict = {}

outname = (".").join(sys.argv[1].split('.')[:-1]) + "_collapsed.fa"

with open(outname, 'w') as outfile:
	for seq in fasta_sequences:
		name, sequence = seq.id, seq.seq.tostring()
		if sequence in seqDict:
			seqDict[sequence].append(name)
		else:
			seqDict[sequence] = [name]
	for sequence in seqDict:
		header = ">" + (" ").join(seqDict[sequence]) + "\n"
		outfile.write(header)
		outfile.write(sequence + "\n")






