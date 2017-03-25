#!/usr/bin/env python3
# Nathan Abell (nsabell)
# GENE211: Genomics
# Final Project

# Load libraries
import sys
import os
import math
import operator
import time
import subprocess
import argparse

# Read in command line arguments with argpase
# Please read the help sections for descriptions of inputs
# The main output is a file called "mpra_test_variants.varInfo.txt" that is fed to generate_oligo_library.py
parser = argparse.ArgumentParser(description='Select Genetic Variants for MPRA Testing from eQTL Summary Statistics')
parser.add_argument('-i', dest='inputFile', required=True, help='Path to summary statistics file')
parser.add_argument('-o', dest='outputFile', required=True, help='Path to output directory (must exist)')
parser.add_argument('-p', dest='pvalue_cutoff', required=False, help='P-value cutoff for significant lead association; default = 10e-7', default = 10**-7)
parser.add_argument('-r', dest='ld_cutoff', required=False, help='R^2 LD Cutoff; default = 1.0 (perfect LD))', default = "1.0")
parser.add_argument('-vcf', dest='vcf', required=False, help='VCF file used for the eQTL study; default = File specific for GENE 211, CEU.chr1-22.n69.p3.MAF01.sorted.vcf', default = "/afs/.ir/users/n/s/nsabell/gene211/ps5/CEU.chr1-22.n69.p3.MAF01.sorted.vcf")
args = parser.parse_args()

# Find the current directory path
main_dir = os.getcwd()

# Main execution
def main():
	'''
	A function to execute the workflow sequentially
	Input: None
	Output: None (just writes to files)
	'''

	# Get all assocation test p-values grouped by gene ID into a dictionary
	associations = load_eqtls(args.inputFile)

	# Find the most strongly associated variants that pass the minimum p-value cutoff
	top_vars = compute_top_variants(associations, args.pvalue_cutoff)

	# Determine if any of the required VCFtools LD calculations are missing
	need_ld_run = False
	for i in range(1,23):
		if not os.path.exists(args.outputFile + "/ld_information/chr" + str(i) + ".list.hap.ld"):
			need_ld_run = True
			break
	
	# Create those files if required based on the tests above
	if need_ld_run:
		unique_top_vars = subset_unique_variants(top_vars)
		generate_ld_file(unique_top_vars, args.outputFile, args.ld_cutoff, args.vcf)

	# Read in the LD R^2 information and store as a dictionary
	ld = load_ld_data(args.outputFile + "/ld_information")

	print(ld[(10,103234839)])

	# Add variants based on whether they pass the LD threshold and write the final set
	filter_by_ld_thresholds(top_vars, ld, args.ld_cutoff, args.outputFile)

	# Recover allele information using the same VCF file to create a varInfo file
	# This calls a short shell script that executes a join command
	subprocess.call(main_dir + "/get_varinfo.sh " + args.vcf +  " " + args.outputFile + "/mpra_test_variants.txt", shell=True)	

def load_eqtls(inputFile):
	'''
	This function parses a set of eQTL summary statistics
	Input: Path to the eQTL data file
	Output: A dictionary where each key is a gene ID and each value is a tuple with the chromosome, position, and -log10 p-value
	'''

	associations = {}

	# Load variant list from file, along with p-values in both sets
	with open(inputFile) as f:
		for line in f:
			data = line.strip().split()

			# Skip header
			if data[0] == "chr":
				continue

			# Add to correct entry, or create if not present
			if data[2] not in associations:
				associations[data[2]] = []

			associations[data[2]].append((int(data[0]), int(data[1]), -1*math.log10(float(data[3]))))

	return(associations)

def compute_top_variants(associations, pvalue_cutoff):
	'''
	This function finds the most associated variants for each gene that pass a p-value threshold
	Input: The output of the load_eqtls function, and a p-value cutoff
	Output: A set containing tuples, each of which has three fields: chromosome, position, and gene ID. 
	These are the most associated variants with each regulated gene.
	'''

	# Create a set container and loop over all eGenes (regulated genes) in the associations dict
	top_vars = set()
	for e in associations:
		
		# Get top variants for each gene
		ranking = sorted(associations[e], key = operator.itemgetter(2), reverse=True)

		# Skip the gene if the best association does not pass the p-value threshold
		if ranking[0][2] < -1*math.log10(pvalue_cutoff):
		 	continue

		 # Add the best variants and any tied variants to the top_vars set
		for cv in ranking:
			if cv[2] == ranking[0][2]:
				top_vars.add((cv[0], cv[1], e))

	return top_vars

def subset_unique_variants(variants):
	'''
	This function just subsets a variant list to only unique positions, even if there are multiple genes
	Input: the output of compute_top_variants
	Returns: another set (as in compute_top_variants) with only unique two-field tuples.
	'''
	vars_uniq = set()
	for var in variants:
		vars_uniq.add((var[0], var[1]))
	return vars_uniq

def generate_ld_file(top_vars, outputFile, ld_cutoff, vcf):
	'''
	This function executes vcftools to compute LD R^2 statistics between all lead "top" variants and those within a 200kb window.
	Input: the top_vars set, the output directory, the ld cutoff, and the vcf file
	Output: a set of written vcftools output files to be parsed later
	'''

	# Write a file containing all the lead variants
	with open(outputFile + "/mpra_lead_variants.txt", "w") as w:
		for tv in top_vars:
			w.write("{0}\t{1}\t{2}\n".format(tv[0], tv[1], int(tv[1])+1))

	# For each of those variants, define a window around it of 100kb in either direction
	with open(outputFile + "/LD_intervals.bed", "w") as w:
		for tv in top_vars:
			w.write("{0}\t{1}\t{2}\n".format(tv[0], int(tv[1])-100000, int(tv[1])+100000))

	# Use vcftools to compute LD of those nearby variants for each population.
	for i in range(1,23):
		subprocess.call("vcftools --vcf {0} --bed {1} --hap-r2-positions {2} --chr {3} --out {4} --min-r2 0.1 &".format(vcf, outputFile + "/LD_intervals.bed", outputFile + "/mpra_lead_variants.txt", str(i), outputFile + "/ld_information/chr" + str(i) + ".list.hap.ld"), shell=True)

def load_ld_data(base_filename):
	'''
	This function loads the LD R^2 information produced by the generate_ld_file function
	Input: the directory containing the vcftools outputs for each chromosome
	Output: a dictionary where each key is a tuple (chromosome, position) and
	each value is a list of tuples (chromosome, position, r^2 value) of all other variants in the window
	'''

	# Initialize the dictionary
	ld = {}

	# Read in each file into the correct directory
	for i in range(1,23):
		with open(base_filename + "/chr" + str(i) + ".list.hap.ld") as f:
			f.readline()
			for line in f:
				data = line.strip().split()
				if "nan" in data[5]:
					continue

				# Add the linkage value to entries for both variants in the pair
				if (int(data[0]), int(data[1])) not in ld:
					ld[(int(data[0]), int(data[1]))] = []
				ld[(int(data[0]), int(data[1]))].append((int(data[2]), int(data[3]), float(data[5])))

				if (int(data[2]), int(data[3])) not in ld:
					ld[(int(data[2]), int(data[3]))] = []
				ld[(int(data[2]), int(data[3]))].append((int(data[0]), int(data[1]), float(data[5])))

	return(ld)

def filter_by_ld_thresholds(top_vars, ld, ld_cutoff, outputFile):
	'''
	This function takes the top lead variants, the LD information, and the cutoffs and produces a final set for synthesis
	Input: the top_vars set, the ld dictionary, the ld cutoff, and the output destination
	Output: a file called "mpra_test_variants.txt" containing the chromosome, position, and 
	'''

	candidate_variants = set()

	for tv in top_vars:

		# Add top variant to our list
		candidate_variants.add((tv[0], tv[1], tv[2]))

		# Add all variants in LD with our top variant to the list
		if (int(tv[0]), int(tv[1])) in ld:
			for variant in ld[(tv[0], tv[1])]:
				if variant[2] >= float(ld_cutoff):

					# Assign any variants added in this way to the associated gene of the linked variant
					candidate_variants.add((variant[0], variant[1], tv[2]))

	# Write everything out to be processed by the shell script
	with open(outputFile + "/mpra_test_variants.txt", "w") as w:
		for cv in candidate_variants:
			w.write("{0}\t{1}\t{2}\n".format(cv[0], cv[1], cv[2]))

if __name__ == "__main__":
	main()


