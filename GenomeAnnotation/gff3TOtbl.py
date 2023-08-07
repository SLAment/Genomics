#!/usr/bin/env python
# encoding: utf-8

# ================== gff3TOtbl =================

# Transform a basic gff3 file to a simple tbl

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

version = 1.0
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
	startcodons_list = ['TTA', 'TTG', 'CTG', 'ATT', 'ATC', 'ATG', 'GTG']
else:
	startcodons_list = ['ATG']


# ---------------------------------
# Read fasta file
# ---------------------------------
seqlens = {}
for seq_record in fastaopen: # Assume there is only one sequence that matters
	mtname = seq_record.id
	seqlens[mtname] = len(seq_record)
	mtseq = seq_record

# ---------------------------------
# Produce the tbl file
# ---------------------------------
def printAttributes(feature):
	for attri in feature.attributes:
		result = ''
		if attri == "ID" or attri == "Name" or attri == 'Parent' or attri == "color":
			pass
		else:
			fullattri = ','.join(feature.attributes[attri])
			result += f"\t\t\t{attri}\t{fullattri}\n"
		sys.stdout.write(result)

for seqid in seqlens.keys():
	genesinctg = [gene for gene in db.features_of_type("gene") if gene.chrom == seqid] # Get only the genes from this contig
	# Head per contig
	print(f">Feature {seqid}")
	print(f"1\t{seqlens[seqid]}\tREFERENCE")

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

		## ---- Check if the first codon is a start codon ----
		partial = False
		cdschildren = [child for child in db.children(gene, featuretype='CDS', order_by='start')]
		if cdschildren != []: # there are CDS features
			if gene.strand == '+':
				firstcodon = cdschildren[0]
				startcodonstart = firstcodon.start - 1 # The -1 is because the gff is base 1
				startcodonseq = mtseq[startcodonstart:startcodonstart+3].seq

				if startcodonseq not in startcodons_list:
					partial = True
					start = "<" + str(start)
			else:
				firstcodon = cdschildren[len(cdschildren) - 1]
				startcodonstart = firstcodon.end - 3
				startcodonseq = mtseq[startcodonstart:startcodonstart+3].seq.reverse_complement()
			
				if startcodonseq not in startcodons_list:
					partial = True
					start = ">" + str(start)
		## ------

		print(f"{start}\t{end}\tgene")
		print(f"\t\t\tlocus_tag\t{gene.id}")
		# print(gene)
		print(f"\t\t\tgene\t{gene.attributes['Name'][0]}")
		printAttributes(gene)

		# Now print the children
		for child in list(db.children(gene, order_by='start', level = 1)): 
			if child.featuretype == "tRNA" or child.featuretype == 'rRNA': # Then print attributes right away
				print(f"{start}\t{end}\t{child.featuretype}")
				printAttributes(child)
		
		if args.mfannot:
			listoffeatures = ['CDS']
		else:
			listoffeatures = ['exon', 'CDS'] # similar to Funannotate's tbl files

		childcount = 0
		for childtype in listoffeatures:
			children = [child for child in db.children(gene, featuretype=childtype, order_by='start')]
			
			if children != []: # tRNA or rRNA might not have children
				if gene.strand == '+':
					# Print the first exon as an mRNA feature
					start = children[0].start
					end = children[0].end

					if partial and not childcount:
						start = "<" + str(start)

					if childtype == "exon":
						print(f"{start}\t{end}\tmRNA") 
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

					if partial and not childcount:
						start = ">" + str(start)

					if childtype == "exon":
						print(f"{start}\t{end}\tmRNA") 
					elif childtype == "CDS":
						print(f"{start}\t{end}\tCDS") 

					childindexes = list(reversed(range(0,len(children))))
		
					# The rest of the exons are just the coordinates
					for exon in childindexes[1:]: # Exclude the last exon
						start = exon.start
						end = exon.end
						print(f"{end}\t{start}") 

			#  --- Print attributes ---
			for mrna in db.children(gene, featuretype='mRNA', order_by='start'): # Assume there is a sinlge mRNA 
				if childtype == "exon":
					for attri in mrna.attributes:
						if attri == "product" or attri == "transcript_id" or attri == "protein_id":
							print(f"\t\t\t{attri}\t{mrna.attributes[attri][0]}")
				elif childtype == "CDS":
					printAttributes(mrna)
			childcount += 1

		# Print the exons and introns too
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
				


