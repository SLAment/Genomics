#!/usr/bin/env python
# encoding: utf-8

# ================== vcfNAreport: Report on missing data levels in vcf files =================
# ==================================================
# Sandra Lorena Ament Velasquez
# Stelkens Lab, Department of Zoology, Stockholm University, Sweden
# 2022/07/14
# +++++++++++++++++++++++++++++++++++++++++++++++++

import argparse # For the fancy options
import sys  # To exit the script, and to pipe out
import os  # For the input name
import re

# ------------------------------------------------------
version = 0.0
versiondisplay = "{0:.2f}".format(version)

# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Report on missing data levels in vcf files *", epilog="") # Create the object using class argparse

# Add options
parser.add_argument('vcf', help="Standard vcf file")
parser.add_argument('--compressed', '-z', help="VCF file is compressed (e.g. with bgzip)", default=False, action='store_true')


try:
	# ArgumentParser parses arguments through the parse_args() method. You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	if args.compressed:
		import gzip # to open compressed files
		vcfopen = gzip.open(args.vcf, 'rt')
	else:
		vcfopen = open(args.vcf, 'r')

except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()

# ------------------------------------------------------
# Define functions
# ------------------------------------------------------
missing1 = re.compile("^\.\/\.:")
missing2 = re.compile("^\.\|\.:")
missing3 = re.compile("^\.:")

def isitNA( covstring):
	if missing1.search(covstring) or missing2.search(covstring) or missing3.search(covstring):
		return True
	else:
		return False

# ------------------------------------------------------
# Read vcf file
# ------------------------------------------------------
basicreport = 0
countofsites = 0
for line in vcfopen:
	if "##" in line: # header
		# sys.stdout.write(line) 
		pass	
	elif "#CHROM" in line: # column names
		# Get samples 
		columns = line.rstrip("\n").split('\t')
		samplesdic = dict( [key, basicreport] for key in columns[9:] )
		samples = list(samplesdic.keys())

	else:
		countofsites += 1
		cols = line.rstrip("\n").split('\t')
		for i in range(0, len(samples)):
			# print(cols[9+i], isitNA(cols[9+i]))
			samplesdic[samples[i]] += isitNA(cols[9+i])

# ------------------------------------------------------
# Report
# ------------------------------------------------------
# Header
print("Sample\tTotal_sites\tMissing_sites\tMissing%")
for sample in samples:
	print(f"{sample}\t{countofsites}\t{samplesdic[sample]}\t{samplesdic[sample]/countofsites:.5f}")

