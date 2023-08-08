#!/usr/bin/env python
# encoding: utf-8

# ================== gff3TOtbl =================

# Transform a basic gff3 file to a simple tbl

# NOTE: The script can't deal with Ontology terms!
# ==================================================
# Sandra Lorena Ament Velasquez
# 2023/07/13
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import os # For the input name
import sys # For reading the input
import argparse # For the fancy options
import gffutils
from Bio import SeqIO
# ------------------------------------------------------

version = 1.21
versiondisplay = "{0:.2f}".format(version)

# ============================
# Check input file
# ============================
# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Transform a gff3 file into a tbl file *", epilog = "It depends on libraries gffutils and biopython.") # Create the object using class argparse

# Add options
parser.add_argument('gff', help="gff3 file")
parser.add_argument('fasta', help="Fasta sequence of the mitochondrial genome")
parser.add_argument('--mfannot', '-m', help="Follow the style of the MFannot tbl file regarding exons and intron", default=False, action='store_true')
parser.add_argument("--code", "-c", help="Genetic code used for traslation (--proteinon) as an NCBI number. Default: 4", default=4, type=int)

# extras
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	fastaopen = SeqIO.parse(args.fasta, "fasta")
	# gffopen = open(args.gff, 'r')
except IOError as msg:  # Check that the file exists
    parser.error(str(msg)) 
    parser.print_help()
# ============================

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
	"rRNA": ["ID", "Name"], 
	"tRNA": ["ID", "Name"]} 

db = gffutils.create_db(data = args.gff, 
	keep_order = True,
	dbfn = dbfnchoice,
	# force = True, # force=True overwrite any existing databases.
	id_spec = id_spec, 
	verbose = False,
	merge_strategy = "create_unique") # Add an underscore an integer at the end for each consecutive occurrence of the same ID 
# print() # if verbose = True

# t1 = time.time()
# from gffutils import inspect; db_results = inspect.inspect(db) # Report
# print()
# print("\n\nIt took {0:.1f}s to create database".format(t1 - t0))
# ---------------------------------

if args.code == 4:
	startcodons_list = ['TTA', 'TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG']
else:
	startcodons_list = ['ATG']

stopcodons_list = ['TAG', 'TAA', 'TGA']

# ---------------------------------
# Read fasta file
# ---------------------------------
records_dict = SeqIO.to_dict(fastaopen)

# ---------------------------------
# Produce the tbl file
# ---------------------------------

def printAttributes(feature): # like funannotate output
	ignoredattributes = ["ID", "Name", 'Parent', "color", "transcript_id", "protein_id"]	
	for attri in feature.attributes:
		if attri in ignoredattributes:
			pass
		elif attri == "product":
			sys.stdout.write(f"\t\t\tproduct\t{feature.attributes['product'][0]}\n")
		else:
			if args.mfannot:
				result = ''
				fullattri = ','.join(feature.attributes[attri])
				result += f"\t\t\t{attri}\t{fullattri}\n"
				sys.stdout.write(result)
			else:
				for item in feature.attributes[attri]:
					if attri == 'Dbxref':
						sys.stdout.write(f"\t\t\tdb_xref\t{item}\n")
					else: 
						sys.stdout.write(f"\t\t\t{attri}\t{item}\n")
	
	if feature.featuretype not in ['gene', 'intron', 'tRNA' , 'rRNA']:
		transcript_id = f"gnl|ncbi|{feature.attributes['ID'][0]}"
		sys.stdout.write(f"\t\t\ttranscript_id\t{transcript_id}\n")
		protein_id = f"gnl|ncbi|{feature.attributes['Parent'][0]}"
		sys.stdout.write(f"\t\t\tprotein_id\t{protein_id}\n")

for seqid in records_dict.keys():
	genesinctg = [gene for gene in db.features_of_type("gene") if gene.chrom == seqid] # Get only the genes from this contig
	ctgseq = records_dict[seqid]

	# Head per contig
	print(f">Feature {seqid}")
	print(f"1\t{len(ctgseq)}\tREFERENCE")

	# Print all genes in this contig
	for gene in genesinctg:
		# In tbl the sense is given by the start and end positions
		if gene.strand == "+":
			start = gene.start
			end = gene.end
		elif gene.strand == "-":
			start = gene.end
			end = gene.start
		else:
			print("The gene {gene.id} has no strand information!")
			sys.exit(1)

		## ---- Make partial genes ----
		partialstart = False
		partialend = False

		cdschildren = [child for child in db.children(gene, featuretype='CDS', order_by='start')]
		if cdschildren != []: # there are CDS features
			if gene.strand == '+':
				## ---- Check if the first codon is a start codon ----
				firstcds = cdschildren[0]
				startcodonstart = firstcds.start - 1 # The -1 is because the gff is base 1
				startcodonseq = ctgseq[startcodonstart:startcodonstart+3].seq

				if startcodonseq not in startcodons_list:
					partialstart = True
					start = "<" + str(start)

				## ---- Check if the last codon is a stop codon ----
				lastcds = cdschildren[len(cdschildren) - 1]
				stopcodonstart = lastcds.end - 3
				stopcodonseq = ctgseq[stopcodonstart:stopcodonstart+3].seq
				
				if stopcodonseq not in stopcodons_list:
					partialend = True
					end = ">" + str(end)

			else:
				## ---- Check if the last codon is a start codon ----
				firstcds = cdschildren[len(cdschildren) - 1]
				startcodonstart = firstcds.end - 3
				startcodonseq = ctgseq[startcodonstart:startcodonstart+3].seq.reverse_complement()
				
				if startcodonseq not in startcodons_list:
					partialstart = True
					start = ">" + str(start)

				## ---- Check if the first codon is a stop codon ----
				lastcds = cdschildren[0]
				stopcodonstart = lastcds.start - 1 # The -1 is because the gff is base 1
				stopcodonseq = ctgseq[stopcodonstart:stopcodonstart+3].seq.reverse_complement()
				
				if stopcodonseq not in stopcodons_list:
					partialend = True
					end = "<" + str(end)

		## ------
		## Gene body
		print(f"{start}\t{end}\tgene")
		if 'Name' in gene.attributes: # often genes have no name
			if gene.attributes['Name'][0] != gene.attributes['ID'][0]: # ignore if they are the same
				print(f"\t\t\tgene\t{gene.attributes['Name'][0]}")
		print(f"\t\t\tlocus_tag\t{gene.id}")
		printAttributes(gene)

		# Now print the children
		for child in list(db.children(gene, order_by='start', level = 1)): # determine the main type of children
			if child.featuretype == "tRNA": # Then print attributes right away because there should be a single tRNA feature
				print(f"{start}\t{end}\t{child.featuretype}")
				printAttributes(child)
			elif child.featuretype == "rRNA":
				if list(db.children(gene, order_by='start', level = 2)) == []: # there are no introns within this rRNA
					print(f"{start}\t{end}\t{child.featuretype}")
			mainchild = child

		listoffeatures = ['exon', 'CDS'] # similar to Funannotate's tbl files (remove 'exon' to have no 'mRNA' feature)

		for childtype in listoffeatures:
			children = [child for child in db.children(gene, featuretype=childtype, order_by='start')]

			if children != []: # tRNA or rRNA might not have children
				if gene.strand == '+':
					# Print the first exon as an mRNA feature
					start = children[0].start
					end = children[0].end

					if partialstart:
						start = "<" + str(start)
					if partialend:
						end = ">" + str(end)

					if childtype == "exon":
						print(f"{start}\t{end}\t{mainchild.featuretype}") 
					elif childtype == "CDS":
						print(f"{start}\t{end}\tCDS")

					# The rest of the exons are just the coordinates
					for exon in children[1:]:
						start = exon.start
						end = exon.end
						print(f"{start}\t{end}") 

				elif gene.strand == '-':
					lastchild = len(children) - 1
					# Print the first exon as an mRNA feature
					start = children[lastchild].end
					end = children[lastchild].start

					if partialstart:
						start = ">" + str(start)

					if childtype == "exon":
						print(f"{start}\t{end}\t{mainchild.featuretype}") 
					elif childtype == "CDS":
						print(f"{start}\t{end}\tCDS") 

					# The rest of the exons are just the coordinates
					childindexes = list(reversed(range(0,len(children))))
					for ind in childindexes[1:]: # Exclude the last exon
						exon = children[ind]
						start = exon.start
						end = exon.end
						print(f"{end}\t{start}") 

			#  --- Print attributes ---
			if childtype == "exon":
				if mainchild.featuretype != 'tRNA': # tRNA was printed above	
					if "protein_id" not in child.attributes or "transcript_id" not in child.attributes:	
						if "product" in child.attributes:
							print(f"\t\t\tproduct\t{mainchild.attributes['product'][0]}")
						
						if mainchild.featuretype != 'rRNA':
							protein_id = f"gnl|ncbi|{mainchild.attributes['Parent'][0]}"
							transcript_id = f"gnl|ncbi|{mainchild.attributes['ID'][0]}"
							sys.stdout.write(f"\t\t\ttranscript_id\t{transcript_id}\n")
							sys.stdout.write(f"\t\t\tprotein_id\t{protein_id}\n")
					else:
						for attri in child.attributes:
							if attri in ["transcript_id", "protein_id"] and mainchild.featuretype != 'rRNA':
								print(f"\t\t\t{attri}\tgnl|ncbi|{child.attributes[attri][0]}")
							elif attri == "product":
								print(f"\t\t\t{attri}\t{child.attributes[attri][0]}")

			elif childtype == "CDS":
				if mainchild.featuretype == 'mRNA': # not necessary for a rRNA
					printAttributes(child)

		# Print the exons and introns too if they exist
		if args.mfannot:
			# Exons
			exonlist = [child for child in db.children(gene, featuretype='exon', order_by='start') if 'orf' not in gene.attributes['Name'][0]]
			count = 0
			# if gene.id == "QC761_0114485": print("hola", [child for child in db.children(gene, order_by='start')])
			for exon in exonlist:
				count += 1
				if gene.strand == '+':
					start = exon.start
					end = exon.end
				else:
					start = exon.end
					end = exon.start					
				print(f"{exon.start}\t{exon.end}\texon")
				print(f"\t\t\tnumber\t{count}")

			# Introns
			intronlist = [child for child in db.children(gene, featuretype='intron', order_by='start')]
			for intron in intronlist:
				if gene.strand == '+':
					start = intron.start
					end = intron.end
				else:
					start = intron.end
					end = intron.start					
				print(f"{intron.start}\t{intron.end}\tintron")
				printAttributes(intron)
				


