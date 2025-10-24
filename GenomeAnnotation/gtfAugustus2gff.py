#!/usr/bin/env python3
# encoding: utf-8

# ================== gtfAugustus2gff =================
# Script to transform the ugly gtf of Augustus into a proper gff3
# ==================================================
# Sandra Lorena Ament Velasquez
# 2025/10/13
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import argparse  # For the fancy options
import os # For the input name
import sys # To exit the script 
import re # Regular expressions
# ------------------------------------------------------
version = 1.1
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* Transform the ugly gtf of Augustus into a proper gff3 *")  # Create the object using class argparse

# Add options
parser.add_argument('GTF', help="GTF file from Augustus")
parser.add_argument("--string", "-s", help="Add a string to the gene model names. Default is to not append anything. Recommended: 'Chr1', but avoid periods (.)!!  ", default='')

# # extras
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	GFFopen = open(args.GTF, 'r')
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()

# ------------------------------------------------------
## Functions

# ------------------------------------------------------
# Dictionary to give individual IDs to CDS/exons in each gene
genedic = {}

sys.stdout.write('##gff-version 3\n')
with open(args.GTF, 'r') as gtffile:
	for line in gtffile:
		contig, source, featuretype, start, end, score, strand, frame, attributes = line.rstrip("\n").split("\t")		# break the line into columns defined by the tab
		
		# Augustus creates many non-standard features, ignore those
		if featuretype == "gene":
			geneid = args.string + attributes
			genedic[geneid] = 1
			attributes = f"ID={geneid};"
			sys.stdout.write(f"{contig}\t{source}\t{featuretype}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attributes}\n")

		elif featuretype == "transcript":
			featuretype = "mRNA"
			parentid = attributes.split('.')[0]
			attributes = f'ID={args.string}{attributes};Parent={args.string}{parentid};'
			sys.stdout.write(f"{contig}\t{source}\t{featuretype}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attributes}\n")

		elif featuretype == "CDS":
			# Make an exon feature
			parentid = args.string + attributes.split(';')[0].split('"')[1]
			geneid = parentid.split('.')[0]
			exonid =  f'{geneid}.exon{genedic[geneid]}'
			attributes = f'ID={exonid};Parent={parentid};'
			
			sys.stdout.write(f"{contig}\t{source}\texon\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attributes}\n")

			# The CDS itself
			cdsid =  f'{geneid}.cds{genedic[geneid]}'
			attributes = f'ID={cdsid};Parent={parentid};'
			sys.stdout.write(f"{contig}\t{source}\tCDS\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attributes}\n")

			# Add this to the count
			genedic[geneid] += 1

