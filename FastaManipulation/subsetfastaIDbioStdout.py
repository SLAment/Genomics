#!/usr/bin/env python
# encoding: utf-8

# ================== Subset fastas using IDs within sequences names =================
# Script to produce a subset of an input fasta file using one or several
# strings that identify specific sequences. 

# It is meant as a faster alternative to subsetfastaID.py but it relies on Biopython.

# But check subsetfasta.py from Lineus work.

# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2014/05/21
# +++++++++++++++++++++++++++++++++++++++++++++++++
# python subsetfastaID.py file.fasta cosito 954

# ------------------------------------------------------
# import argparse # For the fancy options
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import sys # To exit the script 
import os # For the input name
# ------------------------------------------------------
version = 1
versiondisplay = "{0:.2f}".format(version)

# ---------------------------------
# Input from console
# ---------------------------------
try:
	fastafile = sys.argv[1]
	fastaopen = open(fastafile, 'r')
	strings = sys.argv[2:]
except:
	print("Usage: python " + sys.argv[0] + " file.fasta string [ string2 ] ...")
	print("Version" + versiondisplay)
	sys.exit(1)

if strings == []:
	print("Oups! You need to provide one or more strings to subset the fasta file.")
	sys.exit(1)
# ---------------------------------

# ---------------------------------
# Names of sequences in file
# ---------------------------------
# This stores in memory
records_dict = SeqIO.to_dict(SeqIO.parse(fastaopen, "fasta", generic_dna))

for string in strings:
	# Get all sequences that match that string
	matching = [seq for seq in records_dict.keys() if string == seq]

	# Print the fasta sequences of those records
	if len(matching) != 0:
		for seq in matching:
			SeqIO.write(records_dict[seq], sys.stdout, "fasta")
