#!/usr/bin/env python
# encoding: utf-8

# ================== tbl2gff =================

# Transform the tbl into a gff3 file so I can see it in IGV

# v 2.0 - Reshaped the parsing to now obtain all the notes of CDS and introns if present; I also added printing of the last gene which was missing before; it can also deal with rDNA correctly now
# ==================================================
# Sandra Lorena Ament Velasquez
# 2023/07/04
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import os # For the input name
import sys # For reading the input
import argparse # For the fancy options
# ------------------------------------------------------

version = 2.0
versiondisplay = "{0:.2f}".format(version)

# ============================
# Check input file
# ============================
# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Transform a tbl file into a simple gff3 file *") # Create the object using class argparse

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

def dic2string(dic):
	dicstring = ';'
	for key in dic.keys():
		dicstring = dicstring + key + "=" + dic[key] + ";"
	return(dicstring.rstrip(";"))

linecount = 0
for line in tblopen:
	linecount +=1
	if line == '\n' or 'REFERENCE' in line:
		continue
	else:
		cols = line.rstrip("\n").split("\t")		# break the line into columns defined by the tab

	if '>Feature' in line:
		seqid = line.rstrip("\n").split(" ")[1]
	elif 'gene\n' in line:
		if not borja:
			if intronpresent: # print the last intron before the gene, assuming introns are the last thing that comes before a gene
				intronline =  f'{seqid}\tLore\tintron\t{startintron}\t{endintron}\t.\t{sense}\t.\tID={geneID}-intron{introncount};Parent={geneID}-T1' + dic2string(intronattributes)
				genelines.append(intronline)

			if geneName == '':
				geneName = geneID
			print(f'\n{seqid}\tLore\tgene\t{start}\t{end}\t.\t{sense}\t.\tID={geneID};Name={geneName}') # Ugly gff

			if not mrnapresent:
				if trnpresent:
					print(f'{seqid}\tLore\ttRNA\t{start}\t{end}\t.\t{sense}\t.\tID={geneID}-tRNA;Name={geneID};Parent={geneID}' + dic2string(moreattributes) + dic2string(cdsattributes)) # Assume tRNAs always have a product line
				elif rnapresent:
					print(f'{seqid}\tLore\trRNA\t{start}\t{end}\t.\t{sense}\t.\tID={geneID}-rRNA;Name={geneID};Parent={geneID}' + dic2string(moreattributes) + dic2string(cdsattributes)) # Assume there is a single mRNA
				else: # there was no mRNA in the tbl file, so make one for this gene
					print(f'{seqid}\tLore\tmRNA\t{start}\t{end}\t.\t{sense}\t.\tID={geneID}-T1;Name={geneID};Parent={geneID}' + dic2string(moreattributes) + dic2string(cdsattributes)) # Assume there is a single mRNA
			else: # mRNA
				print(f'{seqid}\tLore\tmRNA\t{start}\t{end}\t.\t{sense}\t.\tID={geneID}-T1;Name={geneID};Parent={geneID}' + dic2string(moreattributes) + dic2string(cdsattributes)) # Assume there is a single mRNA

			for child in genelines:
				print(child)

		
		# Re-start the gene
		genelines = []
		geneName = ''
		mrnapresent = False
		trnpresent = False
		cdspresent = False # in Funannotate, the gene comes after CDS
		rnapresent = False
		childcount = 1
		RNAattributes = False
		moreattributes = {}
		cdsattributes = {}

		# Intron stuff
		intronpresent = False
		introncount = 0

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
		RNAattributes = True
	elif 'mRNA\n' in line:
		mrnapresent = True
		trnpresent = False
		RNAattributes = True
	elif 'rRNA\n' in line:
		mrnapresent = False
		trnpresent = False
		rnapresent = True
		RNAattributes = True	
	elif RNAattributes:
		moreattributes[cols[3]] = cols[4]
	# elif 'exon\t' in line: # assuming this always comes after the CDS features, as in the MFannot output (in Funannotate, the gene comes after CDS)
	# 	cdspresent = False
	elif 'intron\n' in line:
		# if there was an intron before, print it
		if intronpresent:
			intronline =  f'{seqid}\tLore\tintron\t{startintron}\t{endintron}\t.\t{sense}\t.\tID={geneID}-intron{introncount};Parent={geneID}-T1' + dic2string(intronattributes)
			genelines.append(intronline)
			intronattributes = {}
		else: # This is the first time we see this intron, so reset
			intronpresent = True
			intronattributes = {}

		if int(cols[0]) > int(cols[1]):
			startintron = int(cols[1])
			endintron = int(cols[0])
			sense = '-'
		else:
			startintron = int(cols[0])
			endintron = int(cols[1])
			sense = '+'

		cdspresent = False # just in case
		introncount += 1

	elif intronpresent: # Notes on that intro
		intronattributes[cols[3]] = cols[4]
	elif cdspresent and len(cols) > 3: # information about the CDS
		cdsattributes[cols[3]] = cols[4]
	elif 'exon\n' in line:
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
			cdspresent = True
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
		
		childcount += 1

# Print the last gene 
print(f'\n{seqid}\tLore\tgene\t{start}\t{end}\t.\t{sense}\t.\tID={geneID};Name={geneName}') # Ugly gff

if not mrnapresent:
	if trnpresent:
		print(f'{seqid}\tLore\ttRNA\t{start}\t{end}\t.\t{sense}\t.\tID={geneID}-tRNA1;Name={geneID};Parent={geneID}' + dic2string(moreattributes) + dic2string(cdsattributes)) # Assume tRNAs always have a product line
	else: # there was no mRNA in the tbl file, so make one for this gene
		print(f'{seqid}\tLore\tmRNA\t{start}\t{end}\t.\t{sense}\t.\tID={geneID}-T1;Name={geneID};Parent={geneID}' + dic2string(moreattributes) + dic2string(cdsattributes)) # Assume there is a single mRNA
else: # mRNA
	print(f'{seqid}\tLore\tmRNA\t{start}\t{end}\t.\t{sense}\t.\tID={geneID}-T1;Name={geneID};Parent={geneID}' + dic2string(moreattributes) + dic2string(cdsattributes)) # Assume there is a single mRNA

for child in genelines:
	print(child)





