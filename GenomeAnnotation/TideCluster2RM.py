#!/usr/bin/env python
# encoding: utf-8

# ================== TideCluster2RM =================
# Script to process the output of TideCluster to produce a consensus library
# for RepeatMasker

# ==================================================
# Sandra Lorena Ament Velasquez
# 2019/03/27 - 2023/07/18
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import argparse  # For the fancy options
import re
# import sys  # To exit the script
# import os  # For the input name
# ------------------------------------------------------
version = 1.1
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* Extract GFF3 features to fasta file *", epilog="The name in the first column of the GFF3 has to be the same as the sequence in the fasta file.")  # Create the object using class argparse

# Add options
parser.add_argument('tarean', help="The tarean_report.tsv file")
parser.add_argument('--superfams', '-s', help="The trc_superfamilies.csv file", default=None, type=str)
parser.add_argument('--monomers', '-m', help="Do NOT duplicate the TRC monomers (consensus), leave them as monomers.", default=False, action='store_true')

# extras
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	args = parser.parse_args()
	tareanopen = open(args.tarean, 'r')
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()

# ------------------------------------------------------

# https://stackoverflow.com/questions/4836710/is-there-a-built-in-function-for-string-natural-sort
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


tabs = [line.rstrip("\n").split("\t") for line in tareanopen]

seqdic = {} # to store the sequences
nameTRC = ''
thiseq = ''
for line in tabs[1:]: # The first line is a header
	if len(line) > 1: # other lines
		counter = 0
		for item in line:
			item = item.replace('"', '')
			counter += 1
			if counter == 2 and 'TRC' in item:
				nameTRC = item
			elif item.startswith('<pre>'):
				thiseq = item
				if item.endswith('<pre>'): # small sequence in the same line
					seqdic[nameTRC] = thiseq
			elif item.endswith('<pre>'):
				thiseq = thiseq + item
				seqdic[nameTRC] = thiseq
	elif len(line) == 1: # it's just a continuation of a sequence
		item = line[0].replace('"', '')
		thiseq = thiseq + item
	# print("------------", nameTRC, thiseq)


superfamilies = {}
if args.superfams is not None:
	sfamtabs = [line.rstrip("\n").split(",") for line in open(args.superfams, 'r')]
	for sftab in sfamtabs[1:]: # I expect a header here
		superfam = sftab[0].replace('"', '')
		clusterID = sftab[1].replace('"', '')
		superfamilies[clusterID] = [superfam]


sorted_names = natural_sort(seqdic.keys())

# Print as a fasta file
for trc in sorted_names:
	# sequence header
	finalseqname = f">{trc}#Satellite"

	if args.superfams is not None:
		if trc in superfamilies.keys():
			finalseqname = f">sf{superfamilies[trc][0]}__{trc}#Satellite"
	print(finalseqname)

	# sequence
	finalseq = f"{seqdic[trc].replace('<pre>', '')}"
	if args.monomers:
		print(finalseq)
	else:
		print(finalseq*2)

