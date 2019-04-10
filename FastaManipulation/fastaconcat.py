#!/usr/bin/env python
# encoding: utf-8

# ================== fastaconcat =================
# Script to concatenate any number of fasta files, based on position at least for now.
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/04/09
# +++++++++++++++++++++++++++++++++++++++++++++++++
version = 1.0

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import sys
from sys import argv

# ============================
# Check input file
# ============================
try:
    fastafile = sys.argv[1]
except:
    print("The program expects the name of the fasta file.\nFor example: \n$ python " + sys.argv[0] + " file.fasta")
    # exit(1)
    sys.exit(1)
# ============================

allfastas = []

# Read all fastas into lists
for fastafile in argv[1:]:
	thisfasta = [seq_record for seq_record in SeqIO.parse(fastafile, "fasta")]
	allfastas.append(thisfasta)

# Get names from the first file
names = [seq.id for seq in allfastas[0]]

# Assume same number of sequences in all files ### STRONG assumption
numseqs = len(names)

concatseqs = []
for i in range(0, numseqs): # for each line of the final fasta
	concatenated = Seq("", generic_dna) # make empty sequence
	for fasta in allfastas: # For each alignment (list of sequences)
		concatenated += fasta[i] # Concatenate the element corresponding to that line of the final fasta
	# Print concatenated sequence
	print(">" + names[i])
	print(concatenated.seq)




