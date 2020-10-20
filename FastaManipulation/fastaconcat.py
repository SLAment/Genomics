#!/usr/bin/env python
# encoding: utf-8

# ================== fastaconcat =================
# Script to concatenate any number of fasta files, based on position at least for now.

# version 2.1 - "from Bio.Alphabet import generic_dna" was removed from BioPython >1.76. See https://biopython.org/wiki/Alphabet
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/04/09
# +++++++++++++++++++++++++++++++++++++++++++++++++
version = 2.1

import argparse # For the fancy options
from Bio.Seq import Seq
from Bio import SeqIO
import sys
from sys import argv

# ============================
# Check input file
# ============================
# try:
#     fastafile = sys.argv[1]
# except:
#     print("The program expects the name of the fasta file.\nFor example: \n$ python " + sys.argv[0] + " file.fasta")
#     # exit(1)
#     sys.exit(1)

# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Concatenate fasta files based on rank or on name *", epilog=".") # Create the object using class argparse

# Add options
# parser.add_argument('fastafile', help="Multifasta file")
parser.add_argument('fastafile', help="One or several multifasta files", type=argparse.FileType('r'), nargs='+') # https://stackoverflow.com/questions/26727314/multiple-files-for-one-argument-in-argparse-python-2-7
parser.add_argument('--name', '-n', help="Concatenate by name (default is by rank)", default=False, action='store_true')


try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
except IOError as msg:  # Check that the file exists
    parser.error(str(msg)) 
    parser.print_help()
# ============================

allfastas = []

# Read all fastas into lists
for fastafile in args.fastafile:
	thisfasta = [seq_record for seq_record in SeqIO.parse(fastafile, "fasta")]
	allfastas.append(thisfasta)

# Get names from the first file
if args.name:
	masternames = {} # A list of all names in all files
	concatlen = 0
	# Collect all names
	for lista in allfastas:
		for seq in lista:
			if seq.id not in masternames.keys():
				masternames[seq.id] = Seq("")

	for lista in allfastas:
		thislen = len(lista[0].seq)

		for name in masternames.keys():
			thisnameseq = "-"*thislen # assume it's not there

			for seq in lista:
				if name == seq.id: # it is there!
					thisnameseq = seq.seq
					break # you found it, so stop looping

			masternames[name] = masternames[name] + thisnameseq

	# Print concatenated sequence
	for seq in masternames.keys():
		print(">" + seq)
		print(masternames[seq])

else:
	names = [seq.id for seq in allfastas[0]]

	# Assume same number of sequences in all files ### STRONG assumption
	numseqs = len(names)

	concatseqs = []
	for i in range(0, numseqs): # for each line of the final fasta
		concatenated = Seq("") # make empty sequence
		for fasta in allfastas: # For each alignment (list of sequences)
			concatenated += fasta[i] # Concatenate the element corresponding to that line of the final fasta
		# Print concatenated sequence
		print(">" + names[i])
		print(concatenated.seq)
