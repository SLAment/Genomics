#!/usr/bin/env python
# encoding: utf-8

# ================== tbl2bed =================

# Transform the tbl into a gff3 file so I can see it in IGV

# ==================================================
# Sandra Lorena Ament Velasquez
# 2023/07/04
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import os # For the input name
import sys # For reading the input
import argparse # For the fancy options
# ------------------------------------------------------

version = 1.0
versiondisplay = "{0:.2f}".format(version)

# ============================
# Check input file
# ============================
# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Transform the tbl back to a BED file *", epilog="BLAST must be locally installed.") # Create the object using class argparse

# Add options
parser.add_argument('tbl', help="tbl file")

# extras
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	tblopen = open(args.tbl, 'r')
except IOError as msg:  # Check that the file exists
    parser.error(str(msg)) 
    parser.print_help()
# ============================

borja = True
cdspresent = False

print('##gff-version 3')

for line in tblopen:
	if line == '\n':
		continue
	else:
		cols = line.rstrip("\n").split("\t")		# break the line into columns defined by the tab

	if '>Feature' in line:
		seqid = line.rstrip("\n").split(" ")[1]
	elif 'gene\n' in line:
		if not borja:
			if geneName == '':
				geneName = geneID
			print(f'{seqid}\tLore\tgene\t{start}\t{end}\t.\t{sense}\t.\tID={geneID};Name={geneName}') # Ugly gff

			if not mrnapresent:

				if trnpresent:
					print(f'{seqid}\tLore\ttRNA\t{start}\t{end}\t.\t{sense}\t.\tID={geneID}-tRNA1;Name={geneID};Parent={geneID};product={thisproduct}') # Assume tRNAs always have a product line
				elif thisproduct: # no mRNA or tRNA, so make a new mRNA feature with product
					print(f'{seqid}\tLore\tmRNA\t{start}\t{end}\t.\t{sense}\t.\tID={geneID}-T1;Name={geneID};Parent={geneID};product={thisproduct}') # Assume there is a single mRNA
				else: # without product
					print(f'{seqid}\tLore\tmRNA\t{start}\t{end}\t.\t{sense}\t.\tID={geneID}-T1;Name={geneID};Parent={geneID}') # Assume there is a single mRNA

			for child in genelines:
				print(child)

		
		# Re-start the gene
		genelines = []
		geneName = ''
		mrnapresent = False
		trnpresent = False
		cdspresent = False
		childcount = 1
		thisproduct = 0

		if cols[0] > cols[1]:
			start = cols[1]
			end = cols[0]
			sense = '-'
		else:
			start = cols[0]
			end = cols[1]
			sense = '+'

		borja = False # not the begining of the file
	elif 'locus_tag' in line: # and newgene:
		geneID = cols[4]
		# print(f'{seqid}\t{start}\t{end}\t{geneID}') # BED file
	elif 'gene\t' in line: # the gene name line
		geneName = cols[4]
	elif 'tRNA\n' in line:
		mrnapresent = False
		trnpresent = True
	elif 'mRNA\n' in line:
		mrnapresent = True
	elif 'product' in line:
		thisproduct = cols[4]
	elif cdspresent and len(cols) > 3:
		cdspresent = False
	elif 'CDS\n' in line or cdspresent:
		if int(cols[0]) > int(cols[1]):
			startcds = int(cols[1])
			endcds = int(cols[0])
			sense = '-'
		else:
			startcds = int(cols[0])
			endcds = int(cols[1])
			sense = '+'

		exonline = f'{seqid}\tLore\texon\t{startcds}\t{endcds}\t.\t{sense}\t.\tID={geneID}-exon{childcount};Parent={geneID}-T1'

		if not cdspresent: # First time reading the CDS
			cdslen = endcds - startcds + 1
			if sense == "+": # Use the previous length to figure out the sense
				phase = 0
			else:
				phase = 3 - (cdslen % 3) # I need to think
		else: # Not the first CDS
			if sense == "+": # Use the previous length to figure out the sense
				phase = 3 - (cdslen % 3)
			else:
				phase = 0 # I need to think

			# Add the new length
			cdslen = cdslen + (endcds - startcds + 1)

		cdsline = f'{seqid}\tLore\tCDS\t{startcds}\t{endcds}\t.\t{sense}\t{phase}\tID={geneID}-cds{childcount};Parent={geneID}-T1'
		
		genelines.append(exonline)
		genelines.append(cdsline)
		cdspresent = True
		childcount += 1





