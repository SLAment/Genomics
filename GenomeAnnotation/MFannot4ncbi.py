#!/usr/bin/env python
# encoding: utf-8

# ================== MFannot4ncbi =================

# MFannot is a useful annotation pipeline for mitochondria. However, the
# output doesn't exactly comply with my needs.

# ==================================================
# Sandra Lorena Ament Velasquez
# 2023/07/02
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import os # For the input name
import sys # For reading the input
import argparse # For the fancy options
from Bio import SeqIO
# ------------------------------------------------------

version = 1.0
versiondisplay = "{0:.2f}".format(version)

# ============================
# Check input file
# ============================
# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Fix MFannot to conform to NCBI requirements better *", epilog="BLAST must be locally installed.") # Create the object using class argparse

# Add options
parser.add_argument('tbl', help="tbl file from MFannot")
parser.add_argument('fasta', help="Fasta sequence of the mitochondrial genome")
parser.add_argument('--locus_tag', '-l', help="Locus tag from NCBI. Default: FUN", default='FUN')
parser.add_argument('--startnumber', '-n', help="Starting number for gene IDs. Default: 1", default='1', type=int)
parser.add_argument('--step', '-s', help="Interval between gene IDs. NCBI suggests a step of 10 in case future annotation updates find new genes in between existing ones. Default: 1", default='1', type=int)
parser.add_argument('--transl_table', '-t', help="Translation table or genetic code. Default: 4", default='4', type=int)

# extras
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	tblopen = open(args.tbl, 'r')
	fastaopen = SeqIO.parse(args.fasta, "fasta")
except IOError as msg:  # Check that the file exists
    parser.error(str(msg)) 
    parser.print_help()
# ============================

for seq_record in fastaopen: # Assume there is only one sequence that matters
	mtname = seq_record.id
	seqlen = len(seq_record)

genecount = args.startnumber # - 1

tRNAdic = {"Alanine": "Ala",
			"Arginine": "Arg",
			"Asparagine": "Asn",
			"Aspartic acid": "Asp",
			"Asparagine": "Asx",
			"Cysteine": "Cys",
			"Glutamine": "Gln",
			"Glutamic acid": "Glu",
			"Glutamine": "Glx",
			"Glycine": "Gly",
			"Histidine": "His",
			"Isoleucine": "Ile",
			"Leucine": "Leu",
			"Lysine": "Lys",
			"Methionine": "Met",
			"Phenylalanine": "Phe",
			"Proline": "Pro",
			"Serine": "Ser",
			"Threonine": "Thr",
			"Tryptophan": "Trp",
			"Tyrosine": "Tyr",
			"Valine": "Val"}

cdsnow = False
for line in tblopen:
	if ">Feature" in line:
		print(f">Feature {mtname}")
		print(f"1\t{seqlen}\tREFERENCE")
	elif 'gene\n' in line:
		if cdsnow: #there was a gene before with CDS
			print(f"\t\t\ttransl_table\t{args.transl_table}")
			cdsnow = False

		sys.stdout.write(line)
		geneID = f'{args.locus_tag}_{"{0:07d}".format(genecount)}'
		print(f'\t\t\tlocus_tag\t{geneID}')
		genecount = genecount + args.step
	elif 'protein_id' in line:
		print(f'\t\t\ttranscript_id\t{geneID}_mrna')
		sys.stdout.write(line.replace('lcl| ', 'gnl|ncbi|'))
	elif 'product\ttransfer RNA' in line:
		anticodon = line.rstrip("\n").replace('\t\t\tproduct\ttransfer RNA ', '')
		print(f'\t\t\tproduct\ttRNA-{tRNAdic[anticodon]}')
	elif 'CDS\n' in line:
		cdsnow = True
		sys.stdout.write(line)
	# elif 'product' in line and cdsnow: # if the codon_start is not specified, then NCBI assumes a codon_start of 1
	# 	print(f"\t\t\tcodon_start\t1")
	# 	sys.stdout.write(line)
	elif 'exon\n' in line and cdsnow: # Assume the exon is right after the CDS
		print(f"\t\t\ttransl_table\t{args.transl_table}")
		cdsnow = False # don't do it in the gene since I did it now
		sys.stdout.write(line)

	else:
		sys.stdout.write(line)




