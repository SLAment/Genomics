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
import os # For the input name
import sys
import gffutils
import datetime
import time
# ------------------------------------------------------
version = 1.1
versiondisplay = "{0:.2f}".format(version)

# Input from console
try:
	gfffile = sys.argv[1]
	geneids = sys.argv[2].split(",")
except:
	print("Usage: python " + sys.argv[0] + " annotation.gff3 gene1_ID,gene2_ID,gene3_ID... > subset.gff3")
	print("Version " + versiondisplay)
	print("Notice: The gene IDs can also be Names, but beware that they might not be unique in a given gff.")
	sys.exit(1)

	gffopen = open(gfffile, 'r')

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
	"similarity": ["Target"], 
	"pseudogene": ["ID", "Name"], 
	"pseudogenic_transcript": ["ID", "Name"], 
	"expressed_sequence_match": ["ID", "Name"]} 

db = gffutils.create_db(data = gfffile, 
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
print('# Subset of ' + sys.argv[1] + ' extracted with GFFSubset.py v. ' + str(versiondisplay) + ' on ' + str(now))


# Parse the gene IDs:
for strid in geneids:
	try: # Try as IDs:
		gene = db[strid]
		printgene(gene)

	except: # Then try gene Name
		# Sometimes the same gene is there for different features
		for gene in db.features_of_type('gene'): # Loop in the database until you fins it
			genename = gene['Name'][0] # same as gene.attributes['Name'][0]
			if strid == genename: 
				printgene(gene)


