#!/usr/bin/env python
# encoding: utf-8

# ================== Psim =================
# Script to calculate Psim sensu Jorda and Kajava (2009) in the program T-REKS
# ==================================================
# Sandra Lorena Ament Velasquez
# 2024/08/14
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import argparse # For the fancy options
from Bio import AlignIO, SeqIO
from Bio.Align import AlignInfo
import sys # To exit the script 
import os # For the input name

# I get a lot of warnings because I use depricated code
import warnings
from Bio import BiopythonDeprecationWarning
warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)
# ------------------------------------------------------
version = 1.0
versiondisplay = "{0:.2f}".format(version)

# ---------------------------------
# Input from console
# ---------------------------------
# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Script to calculate Psim sensu Jorda and Kajava (2009) in the program T-REKS *", epilog=".") # Create the object using class argparse

# Add options
parser.add_argument('fasta', help="Multifasta file")
parser.add_argument('--threshold', '-t', help="The threshold value that is required to add a particular atom (e.g. 0.5). Default: 0", default=0, type=float)

try:
	args = parser.parse_args()
except IOError as msg:  # If no arguments are given
    parser.error(str(msg)) 
    parser.print_help()

# ---------------------------------
# Functions
# ---------------------------------

# Calculate Hamming distance for each sequence
def hamming_distance(seq1, seq2):
	assert len(seq1) == len(seq2)
	return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

# import math
def Psim(alignment):
	# ---- Depricated version, but it works with this biopython version
	# Generate consensus sequence
	summary_align = AlignInfo.SummaryInfo(alignment)
	consensus = summary_align.gap_consensus(threshold=args.threshold) # the variable sites become X
	
	# # print(summary_align)
	# for aln_seq in alignment:
	# 	print(aln_seq.seq)

	print(f">Consensus\n{consensus}")

	# # --- Updated version, but only works with nucleotides
	# # For some reason I can't import Bio.Align from Bio, so I have to read the alignment differently
	# # Read the sequences from a FASTA file
	# alignment = []
	# with open(output.aln) as handle:
	# 	for record in SeqIO.parse(handle, "clustal"):
	# 		print(record)
	# 		alignment.append(record.seq)

	# # Create a Motif object from the sequences
	# motif = create(alignment)

	# # Get the consensus sequence
	# consensus = motif.consensus
	# print(f"Consensus sequence of {wildcards.gene}:", consensus)
	# ---- 
	
	D = [hamming_distance(str(consensus), str(aln_seq.seq)) for aln_seq in alignment]

	# Calculate P_sim
	m = len(alignment)
	l = len(consensus)
	N = m * l

	P_sim = (N - sum(D)) / N

	# for i in range(len(D)):
	# 	# logInvPsim = math.log(1 - ((N - D[i]) / N))
	# 	print(alignment[i].id, D[i])
	# 	# print( 1 - ((N - D[i]) / N) )


	return(P_sim)


alignment = AlignIO.read(args.fasta, "fasta")
P_sim = Psim(alignment)
print(f"Psim\tNo_seqs")
print(f"{P_sim}\t{len(alignment)}")
