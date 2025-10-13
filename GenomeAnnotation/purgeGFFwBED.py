#!/usr/bin/env python
# encoding: utf-8

# ================== purgeGFFwBED =================
# Remove gene models that overlap with a given BED or gff3 file
# ==================================================
# Sandra Lorena Ament Velasquez
# 2025/10/13
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import sys # To exit the script 
import os # For the input name
import argparse  # For the fancy options
import gffutils
import re
# ------------------------------------------------------

version = 1
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* Script to remove gene models that overlap with a given BED file *")  # Create the object using class argparse

# Add options
parser.add_argument('GFF', help="GFF3 file")
parser.add_argument('BED', help="BED file")
parser.add_argument('--gff', '-g', help="The BED file is another gff3 file, e.g. from RepeatMasker.", default=False, action='store_true')
parser.add_argument("--features", "-f", help="String of comma-separated features to use for purging; e.g., TRC_1,TRC_3,TRC_105. Default: use all features in BED.")

# extras
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	GFFopen = open(args.GFF, 'r')
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()
# ============================

# ---------------------------------
# Make database
# ---------------------------------
# This will parse the file, infer the relationships among the features in the file, and store the features and relationships
# See https://pythonhosted.org/gffutils/autodocs/gffutils.create_db.html

dbfnchoice = ':memory:'

# http://daler.github.io/gffutils/database-ids.html
id_spec={"gene": ["ID", "Name"], 
	"mRNA": ["ID", "transcript_id"], 
	"transcript": ["ID", "transcript_id"], 
	"rRNA": ["ID", "Name"],
	"tRNA": ["ID", "Name"],
	"repeat": ["ID", "Name"], 
	"similarity": ["Target"], 
	"pseudogene": ["ID", "Name"], 
	"pseudogenic_transcript": ["ID", "Name"], 
	"expressed_sequence_match": ["ID", "Name"]} 

db = gffutils.create_db(data = args.GFF, 
	keep_order = True,
	dbfn = dbfnchoice,
	# force = True, # force=True overwrite any existing databases.
	id_spec = id_spec, 
	verbose = False,
	merge_strategy = "create_unique") # Add an underscore an integer at the end for each consecutive occurrence of the same ID 

# ---------------------------------

## ----- Deal with the BED file -----
def extract_name(attributes: str) -> str | None:
	"""
	Extract the 'Name' attribute value from a GFF attributes string.
	Example:
	    'ID=TRC_168.000004;Name=TRC_168;Note=308_354;color="#7a1c1a"' → 'TRC_168'
	Returns None if 'Name' is not present.
	"""
	for field in attributes.split(';'):
		if field.startswith('Name='):
			return field.split('=', 1)[1]
	return None

# Prepare the gff or BED file into a list of regions
regionsToPurge = []
with open(args.BED, 'r') as BED:
	tabs = [line.rstrip("\n").split("\t") for line in BED if '##' not in line]

	for line in tabs:
		if args.gff:
			contig = line[0]
			start = line[3]
			end = line[4]

			# Extract identifiers
			name = extract_name(line[8])

		else:
			contig, start, end, name = line[0:4]

		regionsToPurge.append([contig, start, end, name])

if args.features:
	strings = args.features.split(",")
	selectedregions = []

	for region in regionsToPurge:
		if region[3] in strings:
			selectedregions.append(region)

	regionsToPurge = selectedregions

# Sort the regions so they are in order
regionsToPurge.sort(key=lambda x: (x[0], int(x[1])))

## ----- Filter the target gff file -----
# Get iterator of all genes
genes = [gene for gene in db.features_of_type("gene")]


badgenes = []
for gene in genes:
	geneID = gene['ID'][0]
	feature_contig = gene.seqid
		
	for tab in regionsToPurge: # Loop through the regions and see if the gene falls there
		contig, bed_start, bed_end  = tab[:3]

		if feature_contig == contig: #the right contig
			# We want from the start to the stop codon, including the introns, but excluding the UTRs
			allchildren = [child for child in db.children(gene, featuretype='CDS', order_by='start')]

			if len(allchildren) >= 1: # Some genes might not have CDS (e.g. tRNAs)
				feature_start = allchildren[0].start # Start of first CDS to exclude UTRs
				feature_end = allchildren[len(allchildren) - 1].end # End of the last CDS to exclude UTRs

				# Easy case: there is no overlap between the gene and the region
				if (feature_end < int(bed_start)) or (feature_start > int(bed_end)): # not overlapping features!
					pass
				else: 
					badgenes.append(geneID)

				# So BED regions that are *contained* within the genes will be tolerated
				# Genes that are fully contained or overlapping with a region will be removed

badgenes = list(set(badgenes))

# Start an output gff file
sys.stdout.write('##gff-version 3\n')

for gene in genes:
	if gene.id not in badgenes:
		sys.stdout.write(str(gene) + "\n")
		for child in list(db.children(gene)):
			sys.stdout.write(str(child) + "\n")
		# sys.stdout.write("\n")





