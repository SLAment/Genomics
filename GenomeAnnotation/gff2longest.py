#!/usr/bin/env python
# encoding: utf-8

# ================== gff2longest =================
# Script to extract the longest isoform of each protein-coding gene from a 
# GFF3 file and corresponding genome fasta file.

# The script identifies all isoforms (mRNA/transcript features) for each gene,
# calculates their lengths based on CDS or exon content, and outputs only the
# longest isoform for each gene.

# SOURCES
# https://pythonhosted.org/gffutils/autodocs/gffutils.create_db.html
# http://daler.github.io/gffutils/database-ids.html#merge-strategy
# https://github.com/daler/gffutils/blob/master/gffutils/test/test.py
# https://pythonhosted.org/gffutils/database-schema.html
# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

# CAVEATS:
# - The script assumes that the gene is the higher level in the attributes (i.e. polycistronic transcripts are not compatible)
# - It can handle GFF files where mRNA/transcript is the top-level feature (no parent gene)
# - It can't deal with circular features

## VERSION notes
# version 1.0 - Extracts longest isoform as mRNA or CDS, with protein option
# ==================================================
# 2024/04/30
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import argparse  # For the fancy options
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys  # To exit the script
import os  # For the input name
import gffutils
# ------------------------------------------------------

version = 1.0
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* Extract the longest isoform of each protein-coding gene *", epilog="The name in the first column of the GFF3 has to be the same as the sequence in the fasta file. Output is the joined sequence (exons or CDS) of the longest isoform per gene.")

# Add options
parser.add_argument('fastafile', help="Fasta file with the genomic sequences")
parser.add_argument('GFF', help="GFF3 file")
parser.add_argument("--type", "-t", help="The type of sequence to extract: 'mRNA' for joined exons (may include UTRs), 'CDS' for coding sequence only. Default: CDS", default='CDS', choices=['mRNA', 'CDS', 'cds', 'mrna'])
parser.add_argument('--proteinon', '-p', help="Return protein sequences (translates CDS; if --type mRNA, it will still translate based on CDS)", default=False, action='store_true')
parser.add_argument("--code", "-c", help="Genetic code used for translation (--proteinon) as an NCBI number. Default: 1", default=1, type=int)

parser.add_argument("--onlyids", "-n", help="Only keep the ID of the gene/mRNA in the output header", default=False, action='store_true')
parser.add_argument("--onlynames", "-N", help="Only keep the Name of the gene in the output header", default=False, action='store_true')
parser.add_argument("--mRNAids", "-r", help="Use the mRNA/transcript ID instead of the gene ID in the output", default=False, action='store_true')

parser.add_argument('--specificgene', '-g', help="Extract only this specific gene or list of genes separated by commas and no spaces (give ID or Name of gene)")

# extras
parser.add_argument('--output', '-o', help="Name of output file (default is based on input GFF name)", default=None)
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)
parser.add_argument('--verbose', '-b', help="Print a little bit of extra information", default=False, action='store_true')

try:
	args = parser.parse_args()
	GFFopen = open(args.GFF, 'r')
	fastaopen = open(args.fastafile, 'r')
except IOError as msg:
	parser.error(str(msg)) 
	parser.print_help()

# ------------------------------------------------------

# Normalize type argument
if args.type.lower() == "cds": 
	args.type = "CDS"
elif args.type.lower() == "mrna":
	args.type = "mRNA"

# ---------------------------------
# Make database
# ---------------------------------
# This will parse the file, infer the relationships among the features in the file, and store the features and relationships
# See https://pythonhosted.org/gffutils/autodocs/gffutils.create_db.html

dbfnchoice = ':memory:'

id_spec={"gene": ["ID", "Name"], 
	"mRNA": ["ID", "transcript_id"], 
	"transcript": ["ID", "transcript_id"], 
	"rRNA": ["ID", "Name"],
	"tRNA": ["ID", "Name"],
	"pseudogene": ["ID", "Name"], 
	"pseudogenic_transcript": ["ID", "Name"]} 

db = gffutils.create_db(data = args.GFF, 
	keep_order = True,
	dbfn = dbfnchoice,
	# force = True, # force=True overwrite any existing databases.
	id_spec = id_spec, 
	verbose = False,
	merge_strategy = "create_unique") # Add an underscore an integer at the end for each consecutive occurrence of the same ID 

# ---------------------------
### Read the fasta file
# ---------------------------
# This stores in memory
records_dict = SeqIO.to_dict(SeqIO.parse(fastaopen, "fasta"))

# ---------------------------
# Names for input and output
# ---------------------------
if args.output: # User defined
	output_handle = open(args.output, "w")
else:
	input_name = os.path.splitext(args.GFF)[0]
	input_name = os.path.basename(input_name) # Put the output in the working directory, not in the input dir
	suffix = "_longest_" + args.type
	if args.proteinon:
		suffix += "_prot"
	output_handle = open(input_name + suffix + '.fas', "w")

# ---------------------------
### Define some useful functions
# ---------------------------

def seqnamer(seqID, seqname, seq_record, typeseq):
	""" Rename sequence so it's not the chromosome name """
	if args.onlyids and args.onlynames: # both
		seq_record.id = seqID + "_" + seqname
		seq_record.description = ''
	elif args.onlyids: # Name of the output
		seq_record.id = seqID
		seq_record.description = ''
	elif args.onlynames: # Name of the output
		seq_record.id = seqname
		seq_record.description = ''
	else:
		seq_record.id = seqID + '|' + seqname + '|' + typeseq + '|'
		seq_record.description = ''
	return(seq_record)

def get_mrna_type_from_gene(gene, db):
	"""Determine whether children of a gene are 'mRNA' or 'transcript' type"""
	children = db.children(gene, order_by='start', level=1)
	child_types = set()
	for child in children:
		child_types.add(child.featuretype)
	
	if 'mRNA' in child_types: 
		return 'mRNA'
	elif 'transcript' in child_types: 
		return 'transcript'
	return None

def get_isoforms(gene, db, mRNAtype):
	"""Get all mRNA/transcript isoforms for a gene"""
	if mRNAtype:
		return list(db.children(gene, featuretype=mRNAtype, order_by='start'))
	else:
		# If no mRNA/transcript features, return empty list
		return []

def detect_top_level_type(db):
	"""
	Detect whether the GFF has 'gene' as top-level or mRNA/transcript as top-level.
	Returns a tuple: (top_level_type, mrna_type)
	- top_level_type: 'gene', 'mRNA', 'transcript', or None
	- mrna_type: 'mRNA' or 'transcript' (for child features when gene is top-level)
	"""
	# Count features of each type
	gene_count = db.count_features_of_type('gene')
	mrna_count = db.count_features_of_type('mRNA')
	transcript_count = db.count_features_of_type('transcript')
	
	if args.verbose:
		print(f"Feature counts - genes: {gene_count}, mRNAs: {mrna_count}, transcripts: {transcript_count}")
	
	if gene_count > 0:
		# Genes exist, determine mRNA type
		if mrna_count > 0:
			return ('gene', 'mRNA')
		elif transcript_count > 0:
			return ('gene', 'transcript')
		else:
			# Genes exist but no mRNA/transcript children
			return ('gene', None)
	else:
		# No genes, check for top-level mRNA or transcript
		if mrna_count > 0:
			return ('mRNA', 'mRNA')
		elif transcript_count > 0:
			return ('transcript', 'transcript')
		else:
			return (None, None)

def calculate_isoform_length(isoform, db, feature_type='CDS'):
	"""
	Calculate the total length of an isoform based on its CDS or exon features.
	Returns the total length in base pairs.
	"""
	if feature_type == 'CDS':
		children = list(db.children(isoform, featuretype='CDS', order_by='start'))
	else:  # mRNA - use exons
		children = list(db.children(isoform, featuretype='exon', order_by='start'))
	
	total_length = 0
	for child in children:
		total_length += (child.end - child.start + 1)
	
	return total_length

def get_longest_isoform(isoforms, db, feature_type='CDS'):
	"""
	Given a list of isoforms, return the one with the longest CDS or exon content.
	Returns the longest isoform object, or None if no valid isoforms.
	"""
	if len(isoforms) == 0:
		return None
	
	longest = None
	longest_length = 0
	
	for isoform in isoforms:
		length = calculate_isoform_length(isoform, db, feature_type)
		if length > longest_length:
			longest_length = length
			longest = isoform
	
	return longest

def extract_joined_sequence(isoform, db, seq_record, feature_type='CDS'):
	"""
	Extract and join the CDS or exon sequences for an isoform.
	Handles strand orientation (reverse complements minus strand).
	Returns the joined sequence as a Seq object and the strand.
	"""
	if feature_type == 'CDS':
		children = list(db.children(isoform, featuretype='CDS', order_by='start'))
	else:  # mRNA - use exons
		children = list(db.children(isoform, featuretype='exon', order_by='start'))
	
	if len(children) == 0:
		return None, None, None
	
	# Get strand from first child
	strand = children[0].strand
	
	# For CDS, we need to track the phase for proper translation
	first_phase = 0
	last_phase = 0
	
	# Concatenate sequences
	joined_seq = Seq('')
	for i, child in enumerate(children):
		start = child.start - 1  # Convert to 0-based
		stop = child.end
		joined_seq += seq_record.seq[start:stop]
		
		# Track phases for CDS
		if feature_type == 'CDS':
			if i == 0:
				first_phase = int(child.frame) if child.frame != '.' else 0
			last_phase = int(child.frame) if child.frame != '.' else 0
	
	# Handle strand orientation
	if strand == '-':
		# For minus strand, reverse complement
		# Use last_phase to trim the sequence properly (last CDS in genomic order is 5' end of transcript)
		if feature_type == 'CDS' and last_phase > 0:
			joined_seq = joined_seq[:len(joined_seq) - last_phase]
		joined_seq = joined_seq.reverse_complement()
	else:
		# For plus strand, use first_phase to trim
		if feature_type == 'CDS' and first_phase > 0:
			joined_seq = joined_seq[first_phase:]
	
	return joined_seq, strand, first_phase if strand == '+' else last_phase

def extract_cds_for_translation(isoform, db, seq_record):
	"""
	Extract CDS sequence specifically for translation.
	This is used when --type mRNA but --proteinon is set.
	"""
	return extract_joined_sequence(isoform, db, seq_record, feature_type='CDS')

def check_in_focalgenes(feature, featureID, featurename, focalgenes):
	"""Check if a feature (gene or mRNA) matches the focal genes list"""
	if 'Note' in feature.attributes:
		if feature.attributes['Note'][0] in focalgenes:
			return True
	
	if (featureID in focalgenes) or (featurename in focalgenes):
		return True
	
	return False

# ---------------------------
## Only look for some specific gene
# ---------------------------
if args.specificgene: # Only one specific gene, or string of genes, is required
	focalgenes = [gene for gene in args.specificgene.split(",")]

# ---------------------------
## Detect GFF structure
# ---------------------------
top_level_type, global_mrna_type = detect_top_level_type(db)

if args.verbose:
	print(f"Detected structure - Top level: {top_level_type}, mRNA type: {global_mrna_type}")

# ---------------------------
## Process genes and extract longest isoforms
# ---------------------------
gene_count = 0
isoform_count = 0

# Case 1: Gene is the top-level feature
if top_level_type == 'gene':
	for gene in db.features_of_type('gene'):
		geneID = gene['ID'][0]

		## -- Some gffs don't have a name feature (eg. funannotate output)
		try:
			genename = gene['Name'][0]
		except:
			genename = ""
		# --

		# -- Sometimes the first child is transcript and sometimes is mRNA
		# Determine mRNA type for this gene
		mRNAtype = get_mrna_type_from_gene(gene, db)
		# --

		seq_record = records_dict[gene.chrom] # The chromosome sequence
		
		# If the gene is not in the user's list, skip
		if args.specificgene:
			if not check_in_focalgenes(gene, geneID, genename, focalgenes):
				# Also check mRNA IDs if --mRNAids is set
				if args.mRNAids and mRNAtype:
					mrna_children = list(db.children(gene, featuretype=mRNAtype, order_by='start'))
					found = False
					for child in mrna_children:
						if child['ID'][0] in focalgenes:
							found = True
							break
					if not found:
						continue
				else:
					continue

		# Get all isoforms
		isoforms = get_isoforms(gene, db, mRNAtype)
		
		if len(isoforms) == 0:
			# No mRNA/transcript features, try to get CDS/exons directly under gene
			# Treat gene as a single "isoform"
			if args.type == 'CDS':
				cds_children = list(db.children(gene, featuretype='CDS', order_by='start'))
				if len(cds_children) == 0:
					if args.verbose:
						print(f"Gene {geneID} has no CDS features. Skipping.")
					continue
			else:  # mRNA type
				exon_children = list(db.children(gene, featuretype='exon', order_by='start'))
				if len(exon_children) == 0:
					if args.verbose:
						print(f"Gene {geneID} has no exon features. Skipping.")
					continue
			
			# Use gene as the "isoform"
			longest_isoform = gene
			isoformID = geneID
		else:
			# Find the longest isoform
			# Determine which feature type to use for length calculation
			length_type = 'CDS' if args.type == 'CDS' else 'exon'
			longest_isoform = get_longest_isoform(isoforms, db, length_type)
			
			if longest_isoform is None:
				if args.verbose:
					print(f"Gene {geneID} has no valid isoforms with {length_type}. Skipping.")
				continue
			
			isoformID = longest_isoform['ID'][0]
			
			if args.verbose and len(isoforms) > 1:
				print(f"Gene {geneID} has {len(isoforms)} isoforms. Selected {isoformID} as longest.")

		# Extract the sequence
		if args.type == 'CDS':
			joined_seq, strand, phase = extract_joined_sequence(longest_isoform, db, seq_record, 'CDS')
		else:  # mRNA
			joined_seq, strand, phase = extract_joined_sequence(longest_isoform, db, seq_record, 'exon')
		
		if joined_seq is None or len(joined_seq) == 0:
			if args.verbose:
				print(f"Gene {geneID} / isoform {isoformID} has no sequence. Skipping.")
			continue
		
		# Translate to protein if requested
		if args.proteinon:
			if args.type == 'mRNA':
				# Need to get CDS for translation, not exons
				cds_seq, cds_strand, cds_phase = extract_cds_for_translation(longest_isoform, db, seq_record)
				if cds_seq is not None and len(cds_seq) > 0:
					joined_seq = cds_seq.translate(table=args.code)
				else:
					if args.verbose:
						print(f"Gene {geneID} / isoform {isoformID} has no CDS for translation. Skipping.")
					continue
			else:
				joined_seq = joined_seq.translate(table=args.code)
		
		# Create SeqRecord
		output_record = SeqRecord(joined_seq)
		
		# Name the sequence
		if args.mRNAids and isoformID != geneID:
			output_record = seqnamer(isoformID, genename, output_record, args.type)
		else:
			output_record = seqnamer(geneID, genename, output_record, args.type)
		
		# Write to output
		SeqIO.write(output_record, output_handle, "fasta")
		gene_count += 1

# Case 2: mRNA/transcript is the top-level feature (no gene parent)
elif top_level_type in ['mRNA', 'transcript']:
	if args.verbose:
		print(f"Processing {top_level_type} as top-level features (no gene parent found)")
	
	for mrna in db.features_of_type(top_level_type):
		mrnaID = mrna['ID'][0]

		## -- Some gffs don't have a name feature
		try:
			mrnaname = mrna['Name'][0]
		except:
			mrnaname = ""
		# --

		seq_record = records_dict[mrna.chrom] # The chromosome sequence
		
		# If the mRNA is not in the user's list, skip
		if args.specificgene:
			if not check_in_focalgenes(mrna, mrnaID, mrnaname, focalgenes):
				continue

		# Check if this mRNA has CDS (protein-coding)
		cds_children = list(db.children(mrna, featuretype='CDS', order_by='start'))
		if len(cds_children) == 0:
			if args.verbose:
				print(f"mRNA {mrnaID} has no CDS features (not protein-coding). Skipping.")
			continue
		
		# Extract the sequence
		if args.type == 'CDS':
			joined_seq, strand, phase = extract_joined_sequence(mrna, db, seq_record, 'CDS')
		else:  # mRNA
			joined_seq, strand, phase = extract_joined_sequence(mrna, db, seq_record, 'exon')
		
		if joined_seq is None or len(joined_seq) == 0:
			if args.verbose:
				print(f"mRNA {mrnaID} has no sequence. Skipping.")
			continue
		
		# Translate to protein if requested
		if args.proteinon:
			if args.type == 'mRNA':
				# Need to get CDS for translation, not exons
				cds_seq, cds_strand, cds_phase = extract_cds_for_translation(mrna, db, seq_record)
				if cds_seq is not None and len(cds_seq) > 0:
					joined_seq = cds_seq.translate(table=args.code)
				else:
					if args.verbose:
						print(f"mRNA {mrnaID} has no CDS for translation. Skipping.")
					continue
			else:
				joined_seq = joined_seq.translate(table=args.code)
		
		# Create SeqRecord
		output_record = SeqRecord(joined_seq)
		
		# Name the sequence
		output_record = seqnamer(mrnaID, mrnaname, output_record, args.type)
		
		# Write to output
		SeqIO.write(output_record, output_handle, "fasta")
		gene_count += 1

# Case 3: No gene or mRNA features found
else:
	print("ERROR: No gene, mRNA, or transcript features found in the GFF file.")
	print("Available feature types in the database:")
	for featuretype in db.featuretypes():
		count = db.count_features_of_type(featuretype)
		print(f"  {featuretype}: {count}")
	sys.exit(1)

# ---------------------------
## Summary
# ---------------------------
if args.verbose:
	print(f"\nProcessed {gene_count} genes/transcripts.")
	output_name = args.output if args.output else input_name + suffix + '.fas'
	print(f"Output written to: {output_name}")

# Close all opened files
output_handle.close()
fastaopen.close()
GFFopen.close()