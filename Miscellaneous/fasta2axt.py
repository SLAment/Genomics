#!/usr/bin/env python

# ================== fasta2axt =================
# Transform fasta file into an quasi-AXT format for kaks_calculator
# https://github.com/kullrich/kakscalculator2/blob/main/README.md

# I designed the script such that one of the sequences in the fasta file is
# the main focus and all the other sequences are compared against that one.
# ================== fasta2axt =================
# Sandra Lorena Ament Velásquez
# 2025/05/30
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import sys # For reading the input
import os
from Bio import SeqIO # biopython
import argparse # For the fancy options
# ------------------------------------------------------

# ============================
# Make a nice menu for the user
# ============================
version = 1.02
versiondisplay = "{0:.2f}".format(version)

# Make a parser for the options
parser = argparse.ArgumentParser(description="*** fasta2axt ***", epilog="Transform fasta file into an quasi-AXT format for kaks_calculator.") # Create the object using class argparse

# Add basic options
parser.add_argument('fastafile', help="Fasta file containing the query sequences used with MUMmer")
parser.add_argument('-f', '--focus', help="Name of the sequence to be compare all other sequences against. Default is to take the first sequence as focus.", type=str, default = None)

# Useful
parser.add_argument('--ignore', '-i', help="If the focus sequence is not in the alignment, do not throw an error. This is useful if run within Snakemake.", default=False, action='store_true')

# Extras
parser.add_argument('-v', "--version", action='version', version='%(prog)s '+ versiondisplay)


# Parse resulting arguments
args = parser.parse_args()
# ============================

# This stores in memory
records_dict = SeqIO.to_dict(SeqIO.parse(args.fastafile, "fasta"))

focus = ''
if args.focus == None:
	focus = list(records_dict.keys())[0]
else:
	for record in records_dict: # find it in the names of the sequences
		if args.focus == record:
			focus = record
	if focus == '': 
		if args.ignore:
			sys.exit(0)
		else:
			print(f"Sequence {args.focus} not found. Any typos?")
			sys.exit(1)

	records_dict[focus]

focusid = records_dict[focus].id
focusseq = records_dict[focus].seq

for record in records_dict:
	if record == focus:
		pass
	else:
		print(f"{focusid}_vs_{record}")
		print(focusseq)
		print(records_dict[record].seq)
		print()

