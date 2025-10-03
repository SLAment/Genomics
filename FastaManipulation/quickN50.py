#!/usr/bin/env python
# encoding: utf-8

# ================== quickN50 =================
# Calculate a quick N50 value for a given fasta file
# ==================================================
# Sandra Lorena Ament Velasquez
# 2025/10/03
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------

from Bio import SeqIO
import sys

# ------------------------------------------------------
version = 1.1
versiondisplay = "{0:.2f}".format(version)

# ---------------------------------
# Input from console
# ---------------------------------
try:
	fastafile = sys.argv[1]
	# fastaopen = open(fastafile, 'r')
except:
	print("Usage: python quickN50.py file.fasta")
	print("Version" + versiondisplay)
	sys.exit(1)
# ---------------------------------


def calculate_n50(lengths):
    """
    Calculate N50 from a list of contig lengths.
    """
    lengths = sorted(lengths, reverse=True)  # sort from longest to shortest
    total_length = sum(lengths)
    half_length = total_length / 2
    
    running_sum = 0
    for length in lengths:
        running_sum += length
        if running_sum >= half_length:
            return length

def fasta_to_lengths(fasta_file):
    return [len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]


# Produce the N50
contig_lengths = fasta_to_lengths(fastafile)
n50 = calculate_n50(contig_lengths)
print("N50 =", n50)

