#!/usr/bin/env python
# encoding: utf-8

# ================== GFFSubset =================
# Script to extract features of a gff into a new gff based on name or ID. It
# relies on gffutils. For now it assumes that the higher level feature is
# gene.

# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/12/28
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import argparse  # For the fancy options
# import os # For the input name
# import sys
import gffutils
import datetime
import time
# ------------------------------------------------------
version = 2.0
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* Extract features of a gff into a new gff based on name or ID *", epilog="The gene IDs can also be Names, but beware that they might not be unique in a given gff.")  # Create the object using class argparse

# Add options
parser.add_argument('GFF', help="GFF3 file")
parser.add_argument('geneIDs', help="List of genes separated by commas and no spaces (give ID or Name of gene). Eg. gene1_ID,gene2_name,gene3_ID")
parser.add_argument('--keeprepeats', '-r', help="Along with the specified genes, keep also all the repeats in the gff (as from gtfRM2gff.py output)", default=False, action='store_true')

# # extras
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	gffopen = open(args.GFF, 'r')
	geneids = [gene for gene in args.geneIDs.split(",")]
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()

# ------------------------------------------------------

# ---------------------------------
# Make database
# ---------------------------------
# t0 = time.time()
# This will parse the file, infer the relationships among the features in the file, and store the features and relationships
# See https://pythonhosted.org/gffutils/autodocs/gffutils.create_db.html
# I expect a MAKER gff, where CDS have no unique ID

dbfnchoice = ':memory:'

# http://daler.github.io/gffutils/database-ids.html
id_spec={"gene": ["ID", "Name"], 
	"mRNA": ["ID", "transcript_id"], 
	"repeat": ["ID", "Name"], # output of gtfRM2gff.py
	"similarity": ["Target"], 
	"pseudogene": ["ID", "Name"], 
	"pseudogenic_transcript": ["ID", "Name"], 
	"expressed_sequence_match": ["ID", "Name"]} 

db = gffutils.create_db(data = args.GFF, 
	keep_order = True,
	dbfn = dbfnchoice,
	# force = True, # force=True overwrite any existing databases.
	id_spec = id_spec, 
	merge_strategy = "create_unique") # Add an underscore an integer at the end for each consecutive occurrence of the same ID 
	# verbose = True,) 

# t1 = time.time()
# db_results = inspect.inspect(db) # Report
# print("\n\nIt took {0:.1f}s to create database".format(t1 - t0))
# ---------------------------------

## For now let's only work with genes
def printgene(gene):
	print(gene)
	for i in db.children(gene, order_by='start'):
		print(i)

# Print a header
print("##gff-version 3")
now = datetime.datetime.now()
print(f'# Subset of {args.GFF} extracted with GFFSubset.py v. {versiondisplay} on {now}')


# Parse the gene IDs:
for strid in geneids:
	try: # Try as IDs:
		gene = db[strid]
		printgene(gene)

	except: # Then try gene Name
		# Sometimes the same gene is there for different features
		for gene in db.features_of_type('gene'): # Loop in the database until you find it
			genename = gene['Name'][0] # same as gene.attributes['Name'][0]
			if strid in genename: 
				printgene(gene)

# Print in the end also all the repeats in the gff if so desired
if args.keeprepeats:
	for repeat in db.features_of_type('repeat'):
		print(repeat)

