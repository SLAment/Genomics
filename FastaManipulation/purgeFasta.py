#!/usr/bin/env python
# encoding: utf-8

# ================== purgeFasta =================
# Purge a multifasta file from given sequences
# ==================================================
# Sandra Lorena Ament Velasquez
# 2025/03/28
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------

from Bio import SeqIO
import argparse # For the fancy options
import sys

# ------------------------------------------------------
version = 1.1
versiondisplay = "{0:.2f}".format(version)

# ---------------------------------
# Input from console
# ---------------------------------
# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Purge a multifasta file from given sequences *", epilog=".") # Create the object using class argparse

# Add options
parser.add_argument('fasta', help="Multifasta file")
parser.add_argument('list', help="File with the name of the sequences to extract in a column", type=str)
parser.add_argument('--string', '-s', help="The list is a comma-separated string, rather than a one-column file", default=False, action='store_true')
parser.add_argument('--purge', '-v', help="Print the sequences NOT in the list", default=False, action='store_true')

parser.add_argument('--version', "-V", action='version', version='%(prog)s ' + versiondisplay)


try:
	args = parser.parse_args()
except IOError as msg:  # If no arguments are given
    parser.error(str(msg)) 
    parser.print_help()

# ---------------------------------
# Do it!
# ---------------------------------

if args.string:
	tabs = args.list.split(',')
else:
	tabs = {line.rstrip("\n") for line in open(args.list, "r")}  # Read tab file into a set to make it faster to look

for record in SeqIO.parse(args.fasta, "fasta"):
	if args.purge:
		if record.id not in tabs:
			SeqIO.write(record, sys.stdout, "fasta")
	else:
		if record.id in tabs:
			SeqIO.write(record, sys.stdout, "fasta")





