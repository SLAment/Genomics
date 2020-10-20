#!/usr/bin/env python
# encoding: utf-8

# ================== orthogrs2fasta =================
# Script to parse the Orthogroups.csv and SingleCopyOrthogroups.txt output
# files of Orthofinder and to extract each orthogroup in a fasta file
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/05/02
# +++++++++++++++++++++++++++++++++++++++++++++++++

# version 1.1 - "from Bio.Alphabet import generic_dna" was removed from BioPython >1.76. See https://biopython.org/wiki/Alphabet
# ------------------------------------------------------
import sys
import argparse # For the fancy options
from Bio import SeqIO
# ------------------------------------------------------
version = 1.1
versiondisplay = "{0:.2f}".format(version)

# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Parse Orthogroups.csv for Podospora *", epilog="") # Create the object using class argparse

# Add options
parser.add_argument('orthogrscsv', help="The Orthogroups.csv output file of Orthofinder")
parser.add_argument('singlecopies', help="The SingleCopyOrthogroups.txt output file of Orthofinder or equivalent")
parser.add_argument('fastas', help="A space-separated list of fasta files with the original proteins given to OrthoFinder (the names of the files should match the identifier of each sample/sp)", nargs='*') # Take one or more files
parser.add_argument("--outputdir", "-o", help="Path for output directory", default = ".")
parser.add_argument('--samplenames', '-s', help="Replace gene names with the sample names in the output fastas", default=False, action='store_true')

# parser.add_argument("--ref", "-r", help="Reference sample. Default: Podan2", default="Podan2")
# parser.add_argument("--nugrps", "-n", help="Number of orthologs per species per orthogroup. Default: 1", type=int, default=1)
# parser.add_argument("--sample", "-s", help="Sample this number of orthologs. Default: all", type=int, default=float('inf'))

parser.add_argument('--version', "-v", action='version', version='%(prog)s ' + versiondisplay)
parser.add_argument('--verbose', '-b', help="Give some extra information", default=False, action='store_true')

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	orthogrscsvopen = open(args.orthogrscsv, 'r')
	singlecopiesopen = open(args.singlecopies, 'r')
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()

# ------------------------------------------------------
# Parse
# ------------------------------------------------------
tabs = [line.rstrip("\n").split("\t") for line in orthogrscsvopen] 			# Read tab file into a list
onetoones = [line.rstrip("\n") for line in singlecopiesopen] 

samples = tabs[0] # Notice the first field is empty 

orthogroups = {} # Make a master dictionary 

for line in tabs[1:]:
	orthogr = line[0]

	# # Collect names of all proteins in orthogroup
	# if orthogr in onetoones: 
	# 	orthogroups[orthogr] = line[1:]

	# Make a dictionary with nested dics of each species
	if orthogr in onetoones: 
		orthogroups[orthogr] = {} # Make a dictionary for each orthogroup
		for i in range(1,len(samples)): # Ignore the name of the orthogroup
			listgenes = line[i].split(", ")
			if listgenes != ['']: orthogroups[orthogr][samples[i]] = line[i].split(", ") # add nested dictionaries per species (ignore cases without homolog)

# ------------------------------------------------------
# Make fasta 
# ------------------------------------------------------

# Find the right file for each sample
samplefiles = {}
for sample in samples[1:]:
	for fasta in args.fastas:
		if sample + ".fa" in fasta:
			samplefiles[sample] = fasta

# ---- Make a fasta for each orthogroup ----
for ortho in orthogroups:
	# Print some info
	if args.verbose: print("Processing %s ..." % ortho)
	# Open a file for this orthogroup
	ofile = open(args.outputdir + "/" + ortho + ".fa", 'w')

	orthoseqs = [] # Make a list containing the sequences of this orthogroup

	for sample in orthogroups[ortho]:
		# Read the fasta of that sample
		fastaopen = open(samplefiles[sample], 'r')
		records_dict = SeqIO.to_dict(SeqIO.parse(fastaopen, "fasta"))

		# Get the sequence of the right gene
		seq = orthogroups[ortho][sample][0] # Assume there is only one ortholog per species

		# ---------
		if seq not in records_dict.keys(): # The sequence of this sample is not in the keys for some reason
			# Biopython transforms "(" and ")" into "_", so the keys don't match the content on the record_dict
			# Hence, make a new dictionary with those names changed
			newrecords_dict = {}
			for name in records_dict.keys():
				if "(" in name or ")" in name or ":" in name:
					newname = name.replace(')', '_').replace('(', '_').replace(':', '_')
					newrecords_dict[newname] = records_dict[name]
				else:
					newrecords_dict[name] = records_dict[name]

			# Replace old one with new one
			records_dict = newrecords_dict

			# Was that enough?
			if seq not in records_dict.keys():
				print("File %s does not contain the sequence %s" % (samplefiles[sample], seq))
				sys.exit(1)	
		# ---------

		if args.samplenames: # Rename the sequence as the sample
			renamedseq = records_dict[seq]
			renamedseq.id = sample
			renamedseq.description = '' # Removes the description of the protein
			orthoseqs.append(renamedseq)
		else:
			orthoseqs.append(records_dict[seq])

	# Write the ortholog sequences into the new file
	SeqIO.write(orthoseqs, ofile, "fasta")

	# Print some info
	if args.verbose: print("\t... %s completed" % ortho)

