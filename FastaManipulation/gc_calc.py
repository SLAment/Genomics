#!/usr/bin/env python
# encoding: utf-8

# ================== gc_calc =================
# Script to calculate GC content from a (multi) fasta file
# ==================================================
# Sandra Lorena Ament Velasquez
# 2024/10/20
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import argparse # For the fancy options
from Bio import SeqIO

# ------------------------------------------------------
version = 1.0
versiondisplay = "{0:.2f}".format(version)

# ---------------------------------
# Input from console
# ---------------------------------
# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Script to calculate GC content from a (multi) fasta file *", epilog=".") # Create the object using class argparse

# Add options
parser.add_argument('fasta', help="Multifasta file")

try:
	args = parser.parse_args()
except IOError as msg:  # If no arguments are given
    parser.error(str(msg)) 
    parser.print_help()

# ---------------------------------
# Functions
# ---------------------------------

def calculate_gc_content(seq):
	"""
	Calculate GC content for a single sequence, ignoring gaps ('-').
	
	:param seq: The sequence (Bio.Seq object or string) to analyze.
	:return: GC content as a fraction (float).
	"""
	seq_upper = seq.upper()
	gc_count = 0
	valid_bases = 0
	
	for base in seq_upper:
		if base != '-':
			valid_bases += 1
			if base in ['G', 'C']:
				gc_count += 1

	if valid_bases > 0:
		gc_content = gc_count / valid_bases
	else:
		gc_content = 0
	
	return gc_content

print(f"seq\tseq_len\tGC")
for record in SeqIO.parse(args.fasta, "fasta"):
	gc_content = calculate_gc_content(record.seq)
	print(f"{record.id}\t{len(record.seq)}\t{gc_content:.6f}")

