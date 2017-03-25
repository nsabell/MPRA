#!/usr/bin/env python3
# Nathan Abell (nsabell)
# GENE211: Genomics
# Final Project

# Load libraries
import math
import operator
import sys
import os
import time
import pybedtools
import HTSeq as htseq
import subprocess
import argparse
import re
import urllib

# Read in command line arguments with argpase
# Please read the help sections for descriptions of inputs
# The main output is a single file called "MPRA_oligos.fa" and another called "MPRA_oligos_filtered.fa"
parser = argparse.ArgumentParser(description='Generate MPRA-ready oligos from variant information files')
parser.add_argument('-i', dest='inputFile', required=True, help='Path to varInfo file')
parser.add_argument('-w', dest='window', required=False, type=int, default = 150, help='Genomic sequence window length, should be even; default = 150')
parser.add_argument('-r', dest='restrict_seqs', required=False, default = "None", help='Comma-separated list of restriction enzyme sequences to flag ; default = None')
parser.add_argument('-c1', dest='p5_constant', required=False, default = "", help='Constant sequences to add to 5-prime end of each oligo ; default = None')
parser.add_argument('-c2', dest='p3_constant', required=False, default = "", help='Constant sequences to add to 3-prime end of each oligo ; default = None')
parser.add_argument('-o', dest='output_fasta', required=False, default = "MPRA_oligos", help='Output fasta file name ; default = MPRA_oligos.fa')
parser.add_argument('-a', dest='annotation', required=False, default = "gencode.v19.annotation.genes.gff3.gz", help='Annotation to retrieve gene information ; default = Gencode V19 Genes')
parser.add_argument('-n', dest='combnThresh', required=False, type=int, default = 3, help='Maximum number of alternate alleles per oligo ; default = 3')
parser.add_argument('-l', dest='lenThresh', required=False, type=int, default = 160, help='Maximum oligo length ; default = 160')
parser.add_argument('-g', dest='genomePath', required=True, help='Path to genome FASTA from which sequence will be recovered')
args = parser.parse_args()

def main():
	'''
	Main execution
	Input: None, just the command line arguments
	Output: A main FASTA file and a filtered FASTA file (for sequences that are too long of have too many background variants in the same window)
	'''
	
	# Get gene-strand information from gencode v19 (or other)
	strandInfo = readGFF3(args.annotation)

	# Get sequences, metadata, and missing IDs
	sequenceInfo = generateSequences(args.inputFile, int(args.window), strandInfo, args.genomePath,stranded=True)
	writeLibraryFasta(sequenceInfo, args.p5_constant, args.p3_constant, args.output_fasta, args.restrict_seqs, args.combnThresh, args.lenThresh)

def sort_nicely(l):
	'''
	Standard human-friendly sorting function for dictionaries
	Input: a list of alphanumerics
	Output: a sorting that is "human", meaning it increases numerically
	'''
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
	l.sort( key=alphanum_key )
	return(l)


def reverseComplement(string):
	'''
	Standard reverse complement function
	Input: a string
	Output: the reverse complement of that string
	'''

	# Set up the dictionary
	dnaDict = {"A":"T","C":"G","G":"C","T":"A"}
	newString = ""

	# Move backwards to form the RC
	for i in range(0, len(string)):
		newString = dnaDict[string[i].upper()] + newString
	return(newString)

def readGFF3(gff3_path):
	'''
	Function to obtain and read the Gencode v19 GFF3 file using HTSeqs reader function
	Input: the path to the GFF file
	Output: the file is downloaded using urllib, if missing, and read into a HTSeq reader object. Then, the gene IDs and strands are extracted, forming a gene id -> strand dictionary
	'''

	# Download if missing
	if not os.path.exists(gff3_path):
		urllib.request.urlretrieve("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz", gff3_path)

	# Read in with htseq
	gff3 = htseq.GFF_Reader(gff3_path)

	# Form and return the dictionary
	ensg_strand_dict = {}
	for feature in gff3:
		if feature.type == "gene":
			ensg_strand_dict[feature.attr["ID"].split(".")[0]] = feature.iv.strand

	return(ensg_strand_dict)

#Main sequence generating function
#Accounts for haplotype background of arbitrary complexity
def generateSequences(varInfo_path, window, strandDict, genomePath, stranded=False):
	'''
	The main workhorse function: produce the necessary sequences for synthesis
	Input: path to the varInfo file produced with select_mpra_variants.py, the window size, the gene-strand dictionary, the path to the genome fasta, and whether we should care about strand as T/F
	Output: a large dictionary where each key is an ID ("varN") which points to a set of other fields
	'''

	# Start by defining some key parameters, including a counter to keep track of variant ID
	halfwindow = window / 2
	varDict = {}
	counter = 1

	# Keep track of the "last" fields. This is because, for any given a variant, we will read ahead to ensure there are no other variants in the same 150bp window
	lastPosition = 0
	lastStrand = "."
	lastChrom = "chr"

	# Open up the varInfo file
	with open(varInfo_path) as varInfo:

		# Create containers with the other variants in the same window.
		# We add variants to these containers until we get to a variant not in the same window
		# Then, we process what is in the container, and move to the next
		haploVarIDs = []
		haploChroms = []
		haploPositions = []
		haploRefs = []
		haploAlts = []
		haploStrands = []
		haploGenes = []

		# Iterate through the varInfo file
		for line in varInfo:

			# Recover the relevant data
			data = line.split("\t")
			chrom = "chr" + data[0].strip()
			position, gene, ref, alt = int(data[1]), data[2].strip().split(".")[0], data[4].strip(), data[5].strip().split(",")

			# Set strands: either plus strand if stranded=False or if the ID is missing
			# Otherwise, look up in the Gencode annotation-derived dictionary
			if stranded is False:
				strand = "+"
			else: 
				if strandDict.has_key(gene):
					strand = strandDict[gene]
				else:
					strand = "+"

			# If the variant is within the window of the last variant on the same chromosome, or it is the first variant being analyzed, add to the containers
			distance = position - lastPosition + halfwindow
			if 0 < distance < window and chrom == lastChrom:
				haploVarIDs.append(counter)
				haploChroms.append(chrom)
				haploPositions.append(position)
				haploRefs.append(ref)
				haploAlts.append(alt)
				haploStrands.append(strand)
				haploGenes.append([gene])

				# Then update the "last" evaluated variant information
				lastPosition = position
				lastChrom = chrom
				lastStrand = strand
				counter += 1

				# If adding in this way, proceed to the next variant until you get one that is NOT added (because it is too far away)
				continue

			# For each batch of variants that are within 150bp of another variant in this region...
			for i in range(0, len(haploPositions)):

				# Compute the left and right boundaries
				position_left = str(haploPositions[i] - halfwindow)
				position_right = str(haploPositions[i] + halfwindow)

				# Retrieve the sequence information with pybedtools
				interval = " ".join((haploChroms[i], position_left, position_right, ".", ".", "+"))
				sequenceTool = pybedtools.BedTool(interval, from_string=True).sequence(genomePath)
				sequence = sequenceTool.print_sequence().split("\n")[1].upper()

				# Produce the ref and alt oligos, using the information stored in the containers from the varinfo lines
				altOligos = []
				sequence_ref = sequence[0:halfwindow-1] + haploRefs[i] + sequence[halfwindow + len(haploRefs[i]) - 1:]
				for x in range(0,len(haploAlts[i])):
					sequence_alt = sequence[0:halfwindow-1] + haploAlts[i][x] + sequence[halfwindow + len(haploRefs[i]) - 1:]
					altOligos.append([sequence_alt])

				# Make a new entry in the varDict[varID] fields with all the information
				# This is the entry with JUST the centered variant and its alt allele
				# Then, we will take this same sequence, and add in combinatorially all the othe variants in the window
				# Each of these fields shares the same order (i.e. the ref/alt pairs are in the same order)
				varID = "var" + str(haploVarIDs[i])
				varDict[varID] = {}
				varDict[varID]["chrom"] = haploChroms[i]
				varDict[varID]["position"] = haploPositions[i]
				varDict[varID]["reference"] = [sequence_ref]
				varDict[varID]["alternative"] = altOligos
				varDict[varID]["strand"] = haploStrands[i]
				varDict[varID]["haploFlag"] = ["Base"]
				varDict[varID]["refAlleles"] = haploRefs[i]
				varDict[varID]["altAlleles"] = haploAlts[i]
				varDict[varID]["genes"] = ",".join(haploGenes[i])
				varDict[varID]["numOtherVarsInRegion"] = 0

				# Get the positions of all other variants in the same window
				otherPositions = range(0, i) + range(i + 1, len(haploPositions))
				lastOtherPosition = 0

				# Now we go through the other positions in the same window
				for j in otherPositions:

					# Make sure it is within the 150 bp window of this variant
					distance = haploPositions[j] - haploPositions[i] + halfwindow
					if distance == halfwindow or haploPositions[j] == lastOtherPosition:
						continue

					# If it is...
					if 0 < distance < window:

						# Increment the "num other vars in the region" field
						varDict[varID]["numOtherVarsInRegion"] += 1

						# Modify every oligo that has been produced so far
						# These stack up left to right to produe all posible combinations
						for oligoIndex in range(0,len(varDict[varID]["reference"])):

							# Find the position of the new position to modify relative to the center
							# This is necessary for indels
							offset = 0
							offsets = dict(item.split(":") for item in varDict[varID]["haploFlag"][oligoIndex].split(",")[1:]).values()

							for v in offsets:
								pairs = v.split("/")
								offset += len(pairs[1]) - len(pairs[0])
							for allele in haploAlts[j]:
									
								newOligoRef = varDict[varID]["reference"][oligoIndex][0:distance-1+offset] + allele + varDict[varID]["reference"][oligoIndex][distance+offset+len(haploRefs[j])-1:]
								
								varDict[varID]["reference"].append(newOligoRef)

						# Then, for each alternative allele, do the same thing
						for k in range(0,len(varDict[varID]["alternative"])):
							
							mainAlt = varDict[varID]["alternative"][k]
							
							for oligoIndex in range(0,len(mainAlt)):

								offset = 0

								offsets = dict(item.split(":") for item in varDict[varID]["haploFlag"][oligoIndex].split(",")[1:]).values()

								for v in offsets:
									pairs = v.split("/")
									offset += len(pairs[1]) - len(pairs[0])

								if haploPositions[j] - haploPositions[i] > 0:
									offset += len(varDict[varID]["altAlleles"][k]) - len(varDict[varID]["refAlleles"])

								for allele in haploAlts[j]:

									newOligoAlt = mainAlt[oligoIndex][0:distance-1+offset] + allele + mainAlt[oligoIndex][distance+offset+len(haploRefs[j])-1:]
										
									mainAlt.append(newOligoAlt)

									if k == 0:
										varDict[varID]["haploFlag"].append(varDict[varID]["haploFlag"][oligoIndex] + ",var" + str(haploVarIDs[j]) + ":" + haploRefs[j] + "/" + allele)

					lastOtherPosition = haploPositions[j]

			# Finally, since the container is now emptied, the container is reloaded with the variant that that was outside the window of the prior variant
			# Then, the cycle is repeated with each new line
			haploVarIDs = [counter]
			haploChroms = [chrom]
			haploPositions = [position]
			haploRefs = [ref]
			haploAlts = [alt]
			haploStrands = [strand]
			haploGenes = [[gene]]

			lastStrand = strand
			lastPosition = position
			lastChrom = chrom
			counter += 1

	return(varDict)

def restrictionScan(sequence, site):
	'''
	Scan for restriction sites and flag if found
	Allows for N's in the restriction site
	Input: a sequence and a restriction site sequence
	Output: a flag of either "None" or an index of where the idenfied match begins
	'''
	sitePattern = re.compile(site.replace("N","."))
	matches = []
	for match in sitePattern.finditer(sequence):
		matches.append(str(match.start()+1))
	if matches:
		return(",".join(matches))
	return("None")

def writeLibraryFasta(sequenceDict, p5seq, p3seq, out, restSite, combnThresh, lenThresh):
	'''
	Post-process and generate headers
	Write each oligo sequence to standard format fasta file
	Input: the sequence dict produced by generateSequences, the constant 5' and 3' sequences (if they exist),
	the outfile path, restrictions sites, combinatorial thresholds, and length thresholds (on total sequence 
	length)
	'''

	# Open up the output fasta file and the filtered fasta file
	with open(out + ".fa", 'w') as outfile, open(out + "_filtered.fa", 'w') as outfile_filtered:
		
		# Loop through all variants in the sequenceDict, sorted numerically
		for variant in sort_nicely(sequenceDict.keys()):

			# Recover the numbers of background haplotypes and alternative alleles
			varDict = sequenceDict[variant]
			haploNum = len(varDict['haploFlag'])
			altNum = len(varDict['alternative'])

			# Loop through all the background haplotypes
			for i in range(0, haploNum):

				# Get the reference sequence
				referenceSequence = varDict['reference'][i]

				# Reverse complement it if necessary
				if varDict['strand'] == "-":
					referenceSequence = reverseComplement(referenceSequence)

				# If it is too long, trim it down to fit within the provided limit
				if len(referenceSequence) > lenThresh:
					excess = len(referenceSequence) - lenThresh
					trimFromLeft = excess / 2
					trimFromRight = excess - trimFromLeft
					referenceSequence = referenceSequence[trimFromLeft:len(referenceSequence)-trimFromRight]

				# Add the 5' and 3' adapters if necessary
				referenceSequence = p5seq + referenceSequence + p3seq

				# Scan for restriction sites and flag
				if restSite != "None":
					restFlag = restrictionScan(referenceSequence, restSite)
				else:
					restFlag = "NA"

				# If there is a site, make a small mutation in it to salvage it
				if restFlag != "None":
					for restPosition in restFlag.split(","):
						referenceSequence = referenceSequence[0:int(restPosition)+6] + "A" + referenceSequence[int(restPosition)+7:]

				# Format the header with all the fields and metadata separated by underscores
				header = ">" + variant + "_geneIDs=" + varDict['genes'] + "_chrom=" + varDict['chrom'] + "_pos=" + str(varDict['position']) + "_strand=" + varDict['strand'] + "_ref=" + varDict['refAlleles'] + "_alt=" + ','.join(varDict['altAlleles']) + "_allele=ref_haploFlag=" + varDict['haploFlag'][i] + "_restrictionFlag=" + restFlag + "_length=" + str(len(referenceSequence)) + "_numOtherVars=" + str(varDict['numOtherVarsInRegion'])

				# Filter out variants with too maly alternative alleles in the background
				if varDict['numOtherVarsInRegion'] <= combnThresh:
					outfile.write(header + "\n")
					outfile.write(referenceSequence + "\n")				
				else:
					outfile_filtered.write(header + "\n")
					outfile_filtered.write(referenceSequence + "\n")

				# Then, perform the same steps for each alternative allele
				for j in range(0, altNum):
				
					altSequence = varDict['alternative'][j][i]

					# Account for strand
					if varDict['strand'] == "-":
						altSequence = reverseComplement(altSequence)

					# Account for length
					if len(altSequence) > lenThresh:
						excess = len(altSequence) - lenThresh
						trimFromLeft = excess / 2
						trimFromRight = excess - trimFromLeft
						altSequence = altSequence[trimFromLeft:len(altSequence)-trimFromRight]

					# Add adapters
					altSequence = p5seq + altSequence + p3seq

					# Scan for restriction sites
					if restSite != "None":
						restFlag = restrictionScan(altSequence, restSite)
					else:
						restFlag = "None"

					if restFlag != "None":
						for restPosition in restFlag.split(","):
							altSequence = altSequence[0:int(restPosition)+6] + "A" + altSequence[int(restPosition)+7:]

					# Build out header and print
					k = j + 1
					header = ">" + variant + "_geneIDs=" + varDict['genes'] + "_chrom=" + varDict['chrom'] + "_pos=" + str(varDict['position']) + "_strand=" + varDict['strand'] + "_ref=" + varDict['refAlleles'] + "_alt=" + ','.join(varDict['altAlleles']) + "_allele=alt" + str(k) + "_haploFlag=" + varDict['haploFlag'][i] + "_restrictionFlag=" + restFlag + "_length=" + str(len(altSequence)) + "_numOtherVars=" + str(varDict['numOtherVarsInRegion'])

					if varDict['numOtherVarsInRegion'] <= combnThresh:
						outfile.write(header + "\n")
						outfile.write(altSequence + "\n")
					else:
						outfile_filtered.write(header + "\n")
						outfile_filtered.write(altSequence + "\n")

if __name__ == "__main__":
	main()


