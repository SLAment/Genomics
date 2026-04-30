#!/usr/bin/env python
# encoding: utf-8

# ================== gffutils2fasta =================
# Script to extract fasta subsequences out of an input fasta file using a
# corresponding GFF3 file.

# SOURCES
# https://pythonhosted.org/gffutils/autodocs/gffutils.create_db.html
# http://daler.github.io/gffutils/database-ids.html#merge-strategy
# https://github.com/daler/gffutils/blob/master/gffutils/test/test.py
# https://github.com/daler/gffutils/issues/61
# https://pythonhosted.org/gffutils/database-schema.html
# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

# CAVEATS:
# - It can't deal with circular features

## VERSION notes
# version 4.0 - Now handles GFF files where mRNA/transcript is the top-level feature (no parent gene)
# version 3.0 - Now properly handles genes with multiple isoforms (mRNAs/transcripts)
# version 2.33 - Deal with Augustus GFF files that have a transcript rather than mRNA type
# version 2.31 - Now the name of the specific gene will be searched in the Note attribute too
# version 2.3 - New type exoncds gives only the coding part of exons (useful if you have UTRs)
# version 2.2 - Added new option to report all the sequences in the coding sense --insense
# version 2.1 - The database can now take rRNA as a type
# version 2.0 - The database can now take tRNA as a type; the script no longer assumes mRNA is the first level child of gene
# version 1.9 - Now using --onlyids and --onlynames at the same time leads to a combined name
# version 1.7 - Made it more tolerable of CDS without IDs and added an option to extract pseudogenes
# version 1.5 - "from Bio.Alphabet import generic_dna" was removed from BioPython >1.76. See https://biopython.org/wiki/Alphabet
# ==================================================
# Sandra Lorena Ament Velasquez
# 2019/03/27 - 2026/04/30
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import argparse  # For the fancy options
from Bio import SeqIO
from Bio.Seq import Seq
import sys  # To exit the script
import os  # For the input name
import re
import gffutils
# ------------------------------------------------------

version = 4.0
versiondisplay = "{0:.2f}".format(version)
supportedtypes = ["gene", "CDS", "cds", "exon", "exoncds", "noutrs", "similarity", "expressed_sequence_match", "repeat", "pseudogene"] # Unlike the CDS, Exons may contain the UTRs; noutrs is from start to stop codon without introns in nucleotides

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* Extract GFF3 features to fasta file *", epilog="The name in the first column of the GFF3 has to be the same as the sequence in the fasta file.")

# Add options
parser.add_argument('fastafile', help="Fasta file with the genomic sequences")
parser.add_argument('GFF', help="GFF3 file")
parser.add_argument("--type", "-t", help="The feature you want to extract. Options are: "+str(supportedtypes)+", Default: gene. Notice that 'exon' would include the UTRs if present (eg. MAKER output). The option 'noutrs' includes all the body of the gene, except UTRs. The option 'exoncds' gives only the coding part of exons, as opposed to asking for 'CDS', which will trim trailing bases that would shift the frame when translated (if --join is not used).", default='gene', choices = supportedtypes)
parser.add_argument('--proteinon', '-p', help="Return protein sequences (only useful for CDS)", default=False, action='store_true')
parser.add_argument('--join', '-j', help="Join together the features that belong to a single gene (only useful for CDS and exon)", default=False, action='store_true')
parser.add_argument('--specificgene', '-g', help="Extract only this specific gene or list of genes separated by commas and no spaces (give ID or Name of gene)") 
parser.add_argument("--extrabp", "-e", help="Extra base pairs on the side (only works with -t gene, default 0)", type=int, default=0)

parser.add_argument("--onlyids", "-n", help="Only keep the name of the gene ID in the output", default=False, action='store_true')
parser.add_argument("--onlynames", "-N", help="Only keep the name of the gene in the output", default=False, action='store_true')
parser.add_argument("--mRNAids", "-r", help="Use the mRNA ID instead of the gene ID in the output (only for CDS)", default=False, action='store_true')
parser.add_argument("--code", "-c", help="Genetic code used for translation (--proteinon) as an NCBI number. Default: 1", default=1, type=int)

parser.add_argument("--insense", "-s", help="Report all the sequences in the coding sense (experimental!)", default=1, type=int)

# # parser.add_argument("--appendtoname", "-a", help="Append string to the sequence names", default='') # Maybe one day
# # parser.add_argument("--othertype", "-t", help="The feature you want to extract, with your own name", default=NULL) # One day

# extras
parser.add_argument('--output', '-o', help="Name of output file (default is to append 'ed' to the input name)", default=None)
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

# I often write it in small case, so to catch that
if args.type == "cds": 
	args.type = "CDS" # it's necessary in capitals for later
elif args.type == "EXONCDS" or args.type == "exonCDS": 
	args.type = "exoncds"

if args.type == "exoncds" and args.join:
	args.type = "CDS" # It's the same

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
	if args.extrabp > 0:
		output_handle = open(input_name + '_' + str(args.type) + '_p' + str(args.extrabp) +'.fas', "w")
	else:
		output_handle = open(input_name + '_' + str(args.type) + '.fas', "w")

# ---------------------------
### Define some useful functions
# ---------------------------

def getseqbasic(dbobject, seq_record, extrabp = args.extrabp):
	""" Get the sequence of a simple feature """
	start = dbobject.start - 1 - extrabp # The minus one is because in the GFF3 the count starts at 1, but in python it starts at 0s
	stop = dbobject.end + extrabp # Because I use string slice, the last element won't be included
	dbobjectseq = seq_record[start:stop] # Precisely because GFF3 is based 1, so no - 1 is needed
	return(dbobjectseq)

def seqnamer(geneID, genename, geneseq, typeseq = args.type): #, onlyids = args.onlyids, onlynames = args.onlynames):
	""" Rename sequence so it's not the chromosome name """
	if args.onlyids and args.onlynames: # both
		geneseq.id = geneID + "_" + genename
		geneseq.description = ''
	elif args.onlyids: # Name of the output
		geneseq.id = geneID
		geneseq.description = ''
	elif args.onlynames: # Name of the output
		geneseq.id = genename
		geneseq.description = ''
	else:
		geneseq.id = geneID + '|' + genename + '|' + typeseq + '|'	
	return(geneseq)

def get_mrna_type(gene, db):
	"""Determine whether children are 'mRNA' or 'transcript' type"""
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

def check_mrna_in_focalgenes(mrna, mrnaID, mrnaname, focalgenes):
	"""Check if an mRNA matches the focal genes list"""
	if 'Note' in mrna.attributes:
		if mrna.attributes['Note'][0] in focalgenes:
			return True
	
	if (mrnaID in focalgenes) or (mrnaname in focalgenes):
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
## Retrieve the desire features
# ---------------------------
# For pseudogenes
# For pseudogenes
if args.type == "pseudogene": 
	for gene in db.features_of_type('pseudogene'):
		geneID = gene['ID'][0]
		seq_record = records_dict[gene.chrom] # The chromosome sequence
		# Get DNA sequence
		geneseq = getseqbasic(gene, seq_record)
		# Update the naming of the sequence
		geneseq = seqnamer(geneID, geneID, geneseq, typeseq = args.type)

		# Print the sequence 
		SeqIO.write(geneseq, output_handle, "fasta")

### Typical format from a RepeatMasker output
if args.type == "similarity": 
	for gene in db.features_of_type('similarity'):
		geneID = gene['Target'][0]
		seq_record = records_dict[gene.chrom] # The chromosome sequence
		# Get DNA sequence
		geneseq = getseqbasic(gene, seq_record)
		# Update the naming of the sequence
		geneseq = seqnamer(geneID, geneID, geneseq, typeseq = 'repeat')

		# Print the sequence 
		SeqIO.write(geneseq, output_handle, "fasta")

### blastn alignments or the output of `gtfRM2gff.py`
elif (args.type == "expressed_sequence_match") or (args.type == "repeat"):
	for gene in db.features_of_type(args.type):
		geneID = gene['ID'][0]
		genename = gene['Name'][0] # same as gene.attributes['Name'][0]
		seq_record = records_dict[gene.chrom] # The chromosome sequence
		# Get DNA sequence
		geneseq = getseqbasic(gene, seq_record)
		# Update the naming of the sequence
		geneseq = seqnamer(geneID, genename, geneseq, typeseq = args.type)

		# Print the sequence 
		SeqIO.write(geneseq, output_handle, "fasta")

### More canonical genes - now handles both gene-level and mRNA-level top features
else:
	childwarning = False

	# Case 1: Gene is the top-level feature
	if top_level_type == 'gene':
		for gene in db.features_of_type('gene'):
			# dir(gene) --> 'astuple', 'attributes', 'bin', 'calc_bin', 'chrom', 'dialect', 'end', 'extra', 'featuretype', 'file_order', 'frame', 'id', 'keep_order', 'score', 'seqid', 'sequence', 'sort_attribute_values', 'source', 'start', 'stop', 'strand'
			geneID = gene['ID'][0]

			## -- Some gffs don't have a name feature (eg. funannotate output)
			try:
				genename = gene['Name'][0] # same as gene.attributes['Name'][0]
			except:
				genename = "" # An empty string
			# --

			# -- Sometimes the first child is transcript and sometimes is mRNA
			# Determine mRNA type for this gene
			mRNAtype = get_mrna_type(gene, db)
			# --

			seq_record = records_dict[gene.chrom] # The chromosome sequence
			
			# If the gene is not in the user's list, then skip the rest of the code and go to the next gene
			if args.specificgene:
				if 'Note' in gene.attributes:
					if gene.attributes['Note'][0] not in focalgenes: continue

				elif args.mRNAids: # Use the mRNA IDs to search in the focalgenes instead
					mrna_children = list(db.children(gene, featuretype=mRNAtype, order_by='start')) if mRNAtype else []
					found_mrna = False
					for child in mrna_children:
						mrnaID = child['ID'][0]
						if mrnaID in focalgenes:
							found_mrna = True
							break
					if not found_mrna: continue
				elif (geneID not in focalgenes) and (genename not in focalgenes):
					continue

			if args.type == 'gene':
				geneseq = getseqbasic(gene, seq_record) # Get the sequence for this gene
				geneseq = seqnamer(geneID, genename, geneseq, typeseq = args.type) # Update the naming of the sequence

				# Print the sequence 
				SeqIO.write(geneseq, output_handle, "fasta")

			elif args.type == 'noutrs':
				# We want from the start to the stop codon, including the introns, but excluding the UTRs
				# Get all isoforms
				isoforms = get_isoforms(gene, db, mRNAtype)
				
				if len(isoforms) == 0: # Some GFF3 files have no mRNA or transcript features but still have CDS
					# Fallback: get CDS directly under gene
					allchildren = [child for child in db.children(gene, featuretype='CDS', order_by='start')]

					if len(allchildren) >= 1: # Some genes might not even have CDS
						start = allchildren[0].start - 1 - args.extrabp # Start of gene excluding UTRs
						stop = allchildren[len(allchildren) - 1].end + args.extrabp # End of gene excluding UTRs
						geneseq = seq_record[start:stop] # Precisely because GFF3 is based 1, so no - 1 is needed
						geneseq = seqnamer(geneID, genename, geneseq, typeseq = 'gene_noutrs') # Update the naming of the sequence

						# Print the sequence 
						SeqIO.write(geneseq, output_handle, "fasta")
				else: # There is at least one mRNA/transcript
					for isoform in isoforms:
						isoformID = isoform['ID'][0]
						allchildren = [child for child in db.children(isoform, featuretype='CDS', order_by='start')]
						
						if len(allchildren) >= 1: # Some genes might not have CDS
							start = allchildren[0].start - 1 - args.extrabp
							stop = allchildren[len(allchildren) - 1].end + args.extrabp
							geneseq = seq_record[start:stop]
							
							# Use isoform ID in naming if multiple isoforms
							if len(isoforms) > 1 or args.mRNAids:
								geneseq = seqnamer(isoformID, genename, geneseq, typeseq = 'gene_noutrs')
							else:
								geneseq = seqnamer(geneID, genename, geneseq, typeseq = 'gene_noutrs')
							
							SeqIO.write(geneseq, output_handle, "fasta")

			# Getting exons	
			elif args.type == 'exon':
				# Get all isoforms
				isoforms = get_isoforms(gene, db, mRNAtype)
				
				if len(isoforms) == 0:
					# Fallback: get exons directly under gene
					isoforms = [gene]
					fallback_mode = True
				else:
					fallback_mode = False
				
				for isoform in isoforms:
					if fallback_mode:
						isoformID = geneID
						exon_children = db.children(gene, featuretype='exon', order_by='start')
					else:
						isoformID = isoform['ID'][0]
						exon_children = db.children(isoform, featuretype='exon', order_by='start')
					
					child_counter = 1
					child_concat = Seq('')

					for child in exon_children:
						try:
							childID = child['ID'][0]
						except:
							childID = isoformID + '_exon_' + str(child_counter)

						start = child.start - 1
						stop = child.end

						if args.join:
							# Get DNA seq
							child_concat += seq_record[start:stop] # Last letter will be excluded
						else:
							# Get DNA seq
							exonseq = seq_record[start:stop]

							# The naming is a bit different than the normal, so I cannot use the function seqnamer()
							if args.onlyids and args.onlynames:
								exonseq.id = childID + "_" + genename + '_exon.' + str(child_counter)
								exonseq.description = ''
							elif args.onlyids: # Name of the output
								exonseq.id = childID
								exonseq.description = ''
							elif args.onlynames: # Name of the output
								exonseq.id = genename + '_exon.' + str(child_counter)
								exonseq.description = ''
							else:
								exonseq.id = childID + '|' + isoformID + '|exon.' + str(child_counter) + '|'
							
							SeqIO.write(exonseq, output_handle, "fasta")
						
						child_counter += 1

					if args.join and len(child_concat) > 0:
						# Use isoform ID in naming
						if len(isoforms) > 1 or args.mRNAids:
							child_concat = seqnamer(isoformID, genename, child_concat, typeseq = 'concat_exons') # Update the naming of the sequence
						else:
							child_concat = seqnamer(geneID, genename, child_concat, typeseq = 'concat_exons')
						SeqIO.write(child_concat, output_handle, "fasta")

			# Getting CDS	
			elif args.type == 'CDS' or args.type == "exoncds":
				# Get all isoforms
				isoforms = get_isoforms(gene, db, mRNAtype)
				
				if len(isoforms) == 0:
					# Fallback: get CDS directly under gene (some GFFs are structured this way)
					isoforms = [gene]
					fallback_mode = True
				else:
					fallback_mode = False
				
				for isoform in isoforms:
					if fallback_mode:
						isoformID = geneID
						cds_children = list(db.children(gene, featuretype='CDS', order_by='start'))
					else:
						isoformID = isoform['ID'][0]
						cds_children = list(db.children(isoform, featuretype='CDS', order_by='start'))
					
					if len(cds_children) == 0:
						continue  # Skip isoforms without CDS
					
					child_counter = 1
					child_concat = Seq('')
					last_phase = 0  # Track the phase of the last CDS for minus strand

					for child in cds_children:
						# Get the child ID
						if args.mRNAids:
							childID = isoformID
						else:
							# Some gffs don't have an ID for their CDS
							try:
								childID = child['ID'][0]
							except:
								childID = isoformID + '_CDS_' + str(child_counter)
								childwarning = True

						cdsparent = isoformID
						strand = child.strand
						phase = int(child.frame)
						last_phase = phase  # Update for potential use with minus strand

						if (strand != '+') and (strand != '-'):
							print("There is no information of strand direction!")
							sys.exit(1)

						# Produce a protein sequence for the CDS in each exon (taking care of the phases)
						if not args.join:
							if strand == '+':
								start = child.start - 1 + phase
								raw_stop = child.end

								if args.type == "CDS":
									# Trim the extra bases at the end that are lost because the CDS are not joint
									howmanycodons = (raw_stop - start) // 3
									stop = (howmanycodons * 3) + start 
								elif args.type == "exoncds":
									start = child.start - 1
									stop = raw_stop

								# Get DNA seq and translate
								cdsseq = seq_record[start:stop]

								# If protein sequences are required
								if args.proteinon:
									cdsseq.seq = cdsseq.seq.translate(table = args.code)

							elif strand == '-': 
								raw_start = child.start - 1
								stop = child.end - phase 

								if args.type == "CDS":
									# Trim the extra bases at the end that are lost because the CDS are not joint
									howmanycodons = (stop - raw_start) // 3
									start = stop - (howmanycodons * 3)
								elif args.type == "exoncds":
									start = raw_start
									stop = child.end

								if args.proteinon or args.type == "exoncds": # Get DNA seq, reverse complement, and translate
									cdsseq = seq_record[start:stop].reverse_complement()
									cdsseq.id = seq_record.id
									cdsseq.description = seq_record.description
									if args.proteinon:
										cdsseq.seq = cdsseq.seq.translate(table = args.code)
								else:
									cdsseq = seq_record[start:stop]

							if args.onlyids and args.onlynames:
								cdsseq.id = childID + "_" + genename + '_CDS.' + str(child_counter)
								cdsseq.description = ''
							elif args.onlyids:
								cdsseq.id = childID
								cdsseq.description = ''
							elif args.onlynames:
								cdsseq.id = genename + '_CDS.' + str(child_counter)
								cdsseq.description = ''
							else:
								cdsseq.id = childID + '|' + cdsparent + '|CDS.' + str(child_counter) + '|'
						
							SeqIO.write(cdsseq, output_handle, "fasta")
							
							child_counter += 1
						
						# Produce a single CDS for the entire isoform
						else:
							if strand == '+':
								if child_counter == 1: # Consider it only on the beginning of the protein
									start = child.start - 1 + phase
									child_counter += 1
								else: # All the other cds should be in frame
									start = child.start - 1

							elif strand == '-':
								start = child.start - 1
								child_counter += 1 # OJO

							stop = child.end
							# Get DNA seq
							child_concat += seq_record[start:stop]

					# Print joined CDS for this isoform
					if args.join and args.type == 'CDS' and len(child_concat) > 0:
						try: # Sometimes genes won't have a CDS feature, in which case we better move on...
							# Get the strand from the last CDS child processed
							if strand == '+':
								if args.proteinon:
									cdsseq_concat = child_concat.seq.translate(table = args.code)
									child_concat.seq = cdsseq_concat
							elif strand == '-':
								# Use phase from the last CDS (which is the 5' end in minus strand) and reverse complement it
								child_concat = child_concat[:len(child_concat) - last_phase].reverse_complement()
								# Edit Seq object to have the same names again
								child_concat.id = seq_record.id
								child_concat.description = seq_record.description

								if args.proteinon: # Translate to protein
									cdsseq_concat = child_concat.seq.translate(table = args.code)
									child_concat.seq = cdsseq_concat

							elif args.insense and strand == '-': # It's not a protein but make it in the right sense
								child_concat = child_concat[:len(child_concat) - last_phase].reverse_complement()

							# Use isoform ID in naming
							if args.mRNAids or len(isoforms) > 1:
								child_concat = seqnamer(isoformID, genename, child_concat, typeseq = args.type)
							else:
								child_concat = seqnamer(geneID, genename, child_concat, typeseq = args.type)

							SeqIO.write(child_concat, output_handle, "fasta")

						except Exception as e:
							if args.verbose: 
								print(f"The isoform {isoformID} of gene {geneID} {genename} is problematic: {e}. Maybe it doesn't have a CDS? Skipped.")

	# Case 2: mRNA/transcript is the top-level feature (no gene parent)
	elif top_level_type in ['mRNA', 'transcript']:
		if args.verbose:
			print(f"Processing {top_level_type} as top-level features (no gene parent found)")
		
		for mrna in db.features_of_type(top_level_type):
			# dir(mrna) --> 'astuple', 'attributes', 'bin', 'calc_bin', 'chrom', 'dialect', 'end', 'extra', 'featuretype', 'file_order', 'frame', 'id', 'keep_order', 'score', 'seqid', 'sequence', 'sort_attribute_values', 'source', 'start', 'stop', 'strand'
			mrnaID = mrna['ID'][0]

			## -- Some gffs don't have a name feature
			try:
				mrnaname = mrna['Name'][0] # same as mrna.attributes['Name'][0]
			except:
				mrnaname = "" # An empty string
			# --

			seq_record = records_dict[mrna.chrom] # The chromosome sequence
			
			# If the mRNA is not in the user's list, then skip the rest of the code and go to the next mRNA
			if args.specificgene:
				if not check_mrna_in_focalgenes(mrna, mrnaID, mrnaname, focalgenes):
					continue

			if args.type == 'gene':
				# When user asks for 'gene' but only mRNA exists, extract the mRNA span
				geneseq = getseqbasic(mrna, seq_record) # Get the sequence for this mRNA
				geneseq = seqnamer(mrnaID, mrnaname, geneseq, typeseq = 'mRNA') # Update the naming of the sequence

				# Print the sequence 
				SeqIO.write(geneseq, output_handle, "fasta")

			elif args.type == 'noutrs':
				# We want from the start to the stop codon, including the introns, but excluding the UTRs
				allchildren = [child for child in db.children(mrna, featuretype='CDS', order_by='start')]

				if len(allchildren) >= 1: # Some mRNAs might not have CDS
					start = allchildren[0].start - 1 - args.extrabp # Start of mRNA excluding UTRs
					stop = allchildren[len(allchildren) - 1].end + args.extrabp # End of mRNA excluding UTRs
					geneseq = seq_record[start:stop] # Precisely because GFF3 is based 1, so no - 1 is needed
					geneseq = seqnamer(mrnaID, mrnaname, geneseq, typeseq = 'gene_noutrs') # Update the naming of the sequence

					# Print the sequence 
					SeqIO.write(geneseq, output_handle, "fasta")

			# Getting exons	
			elif args.type == 'exon':
				exon_children = db.children(mrna, featuretype='exon', order_by='start')
				
				child_counter = 1
				child_concat = Seq('')

				for child in exon_children:
					try:
						childID = child['ID'][0]
					except:
						childID = mrnaID + '_exon_' + str(child_counter)

					start = child.start - 1
					stop = child.end

					if args.join:
						# Get DNA seq
						child_concat += seq_record[start:stop] # Last letter will be excluded
					else:
						# Get DNA seq
						exonseq = seq_record[start:stop]

						# The naming is a bit different than the normal, so I cannot use the function seqnamer()
						if args.onlyids and args.onlynames:
							exonseq.id = childID + "_" + mrnaname + '_exon.' + str(child_counter)
							exonseq.description = ''
						elif args.onlyids: # Name of the output
							exonseq.id = childID
							exonseq.description = ''
						elif args.onlynames: # Name of the output
							exonseq.id = mrnaname + '_exon.' + str(child_counter)
							exonseq.description = ''
						else:
							exonseq.id = childID + '|' + mrnaID + '|exon.' + str(child_counter) + '|'
						
						SeqIO.write(exonseq, output_handle, "fasta")
					
					child_counter += 1

				if args.join and len(child_concat) > 0:
					child_concat = seqnamer(mrnaID, mrnaname, child_concat, typeseq = 'concat_exons') # Update the naming of the sequence
					SeqIO.write(child_concat, output_handle, "fasta")

			# Getting CDS	
			elif args.type == 'CDS' or args.type == "exoncds":
				cds_children = list(db.children(mrna, featuretype='CDS', order_by='start'))
				
				if len(cds_children) == 0:
					continue  # Skip mRNAs without CDS
				
				child_counter = 1
				child_concat = Seq('')
				last_phase = 0  # Track the phase of the last CDS for minus strand

				for child in cds_children:
					# Get the child ID
					if args.mRNAids:
						childID = mrnaID
					else:
						# Some gffs don't have an ID for their CDS
						try:
							childID = child['ID'][0]
						except:
							childID = mrnaID + '_CDS_' + str(child_counter)
							childwarning = True

					cdsparent = mrnaID
					strand = child.strand
					phase = int(child.frame)
					last_phase = phase  # Update for potential use with minus strand

					if (strand != '+') and (strand != '-'):
						print("There is no information of strand direction!")
						sys.exit(1)

					# Produce a protein sequence for the CDS in each exon (taking care of the phases)
					if not args.join:
						if strand == '+':
							start = child.start - 1 + phase
							raw_stop = child.end

							if args.type == "CDS":
								# Trim the extra bases at the end that are lost because the CDS are not joint
								howmanycodons = (raw_stop - start) // 3
								stop = (howmanycodons * 3) + start 
							elif args.type == "exoncds":
								start = child.start - 1
								stop = raw_stop

							# Get DNA seq and translate
							cdsseq = seq_record[start:stop]

							# If protein sequences are required
							if args.proteinon:
								cdsseq.seq = cdsseq.seq.translate(table = args.code)

						elif strand == '-': 
							raw_start = child.start - 1
							stop = child.end - phase 

							if args.type == "CDS":
								# Trim the extra bases at the end that are lost because the CDS are not joint
								howmanycodons = (stop - raw_start) // 3
								start = stop - (howmanycodons * 3)
							elif args.type == "exoncds":
								start = raw_start
								stop = child.end

							if args.proteinon or args.type == "exoncds": # Get DNA seq, reverse complement, and translate
								cdsseq = seq_record[start:stop].reverse_complement()
								cdsseq.id = seq_record.id
								cdsseq.description = seq_record.description
								if args.proteinon:
									cdsseq.seq = cdsseq.seq.translate(table = args.code)
							else:
								cdsseq = seq_record[start:stop]

						if args.onlyids and args.onlynames:
							cdsseq.id = childID + "_" + mrnaname + '_CDS.' + str(child_counter)
							cdsseq.description = ''
						elif args.onlyids:
							cdsseq.id = childID
							cdsseq.description = ''
						elif args.onlynames:
							cdsseq.id = mrnaname + '_CDS.' + str(child_counter)
							cdsseq.description = ''
						else:
							cdsseq.id = childID + '|' + cdsparent + '|CDS.' + str(child_counter) + '|'
					
						SeqIO.write(cdsseq, output_handle, "fasta")
						
						child_counter += 1
					
					# Produce a single CDS for the entire mRNA
					else:
						if strand == '+':
							if child_counter == 1: # Consider it only on the beginning of the protein
								start = child.start - 1 + phase
								child_counter += 1
							else: # All the other cds should be in frame
								start = child.start - 1

						elif strand == '-':
							start = child.start - 1
							child_counter += 1 # OJO

						stop = child.end
						# Get DNA seq
						child_concat += seq_record[start:stop]

				# Print joined CDS for this mRNA
				if args.join and args.type == 'CDS' and len(child_concat) > 0:
					try: # Sometimes mRNAs won't have a CDS feature, in which case we better move on...
						# Get the strand from the last CDS child processed
						if strand == '+':
							if args.proteinon:
								cdsseq_concat = child_concat.seq.translate(table = args.code)
								child_concat.seq = cdsseq_concat
						elif strand == '-':
							# Use phase from the last CDS (which is the 5' end in minus strand) and reverse complement it
							child_concat = child_concat[:len(child_concat) - last_phase].reverse_complement()
							# Edit Seq object to have the same names again
							child_concat.id = seq_record.id
							child_concat.description = seq_record.description

							if args.proteinon: # Translate to protein
								cdsseq_concat = child_concat.seq.translate(table = args.code)
								child_concat.seq = cdsseq_concat

						elif args.insense and strand == '-': # It's not a protein but make it in the right sense
							child_concat = child_concat[:len(child_concat) - last_phase].reverse_complement()

						child_concat = seqnamer(mrnaID, mrnaname, child_concat, typeseq = args.type)

						SeqIO.write(child_concat, output_handle, "fasta")

					except Exception as e:
						if args.verbose: 
							print(f"The mRNA {mrnaID} {mrnaname} is problematic: {e}. Maybe it doesn't have a CDS? Skipped.")

	# Case 3: No gene or mRNA features found
	else:
		print("WARNING: No gene, mRNA, or transcript features found in the GFF file.")
		print("Available feature types in the database:")
		for featuretype in db.featuretypes():
			count = db.count_features_of_type(featuretype)
			print(f"  {featuretype}: {count}")

	if childwarning and args.onlyids and not args.join: 
		print("WARNING: The CDS have no IDs so the sequences will have no names!")


# Close all opened files
output_handle.close()
fastaopen.close()
GFFopen.close()

"""
gene --> A region (or regions) that includes all of the sequence elements necessary to encode a functional transcript.
exons --> A region of the transcript sequence within a gene which is not removed from the primary RNA transcript by RNA splicing.
CDS --> A contiguous sequence which begins with, and includes, a start codon and ends with, and includes, a stop codon.

A GFF3 file has the following fields

##gff-version 3
0	ref_seq	Thamnoliak127a3n4s400-scf_the-scaffold
1	source	SNAP
2	type	gene	( This is constrained to be either: (a) a term from the "lite" sequence ontology, SOFA; or (b) a SOFA accession number. The latter alternative is distinguished using the syntax SO:000000.)	
3	start_position	835
4	stop_position	1002
5	score	16.007
6	strand	+
7	phase_codon	.
8	attributes	ID=snap.gene.000001; Name=Thamnoliak127a3n4s400-scf_the-scaffold-snap.1
"""

# ---------------------------------
# Learning how to use gffutils 
# ---------------------------------
# What features does it have?
# print(db_results)
# print(list(db_results['featuretype'].keys()))
# allfeats = list(db.all_features()) # Iterator of all features
# allIDsallfeats = [f.id for f in db.all_features()]
# print(allIDsallfeats)

# # Print the whole thing
# for gene in db.features_of_type('gene'):
# 	print(gene)
# 	for child in list(db.children(gene)):
# 		print(child)
# --------------------------------------------

# print("--- %s seconds ---" % (time.time() - start_time))
