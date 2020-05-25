#!/usr/bin/env python
# encoding: utf-8

# ================== nexus2markers =================
# Script to extract all partitions in a simple nexus file as individual fasta files
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020/05/25
# +++++++++++++++++++++++++++++++++++++++++++++++++
# https://biopython.org/wiki/AlignIO
# http://biopython.org/DIST/docs/api/Bio.Align.MultipleSeqAlignment-class.html
# http://biopython.org/DIST/docs/api/Bio.Nexus.Nexus-pysrc.html
# https://gist.github.com/brantfaircloth/2999578
# https://biopython.org/wiki/Concatenate_nexus

version = 1

from Bio.Nexus import Nexus
import sys
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO

# ============================
# Check input file
# ============================
try:
    nexfile = sys.argv[1]
except:
    print("The program expects the name of the nexus file.\nFor example: \n$ python " + sys.argv[0] + " file.nex")
    sys.exit(1)

# ============================

# Read nexus file with Bio.Nexus
master = Nexus.Nexus() # Create a nexus object
master.read(nexfile) # update

markers = list(master.charsets.keys()) # Get list of partitions
# lenmaster = master.nchar

for marker in markers:
	# Get positions included in this partition
	partitioncols = master.charsets[marker]

	# Start an empty msa object
	newalign = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))

	# For each taxon, get those columns:
	for taxon in master.taxlabels:
		record = master.matrix[taxon] # Get sequence of this taxon

		newseq = "" # Start an empty new sequence
		for column in partitioncols:
			newseq += record[column]
		newalign.add_sequence(taxon, newseq) # Add the sliced sequenced to the new msa

	# Write a fasta file for each marker in the working directory
	output_handle = open(marker + ".fas", "w")
	SeqIO.write(newalign, output_handle, "fasta")
