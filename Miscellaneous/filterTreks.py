#!/usr/bin/env python

# ================== filterTreks =================
# Filter the output of the T-reks program (https://bioinfo.crbm.cnrs.fr/index.php?route=tools&tool=3)
# ================== fasta2axt =================
# Sandra Lorena Ament Velásquez
# 2025/06/12-20
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import sys # For reading the input
import os
import argparse # For the fancy options
# ------------------------------------------------------

# ============================
# Make a nice menu for the user
# ============================
version = 1.0
versiondisplay = "{0:.2f}".format(version)

# Make a parser for the options
parser = argparse.ArgumentParser(description="*** filterTreks ***", epilog="Filter the output of the T-reks program.") # Create the object using class argparse

# Add basic options
parser.add_argument('treks', help="Output table of the T-reks program")
parser.add_argument('-m', '--minlen', help="Minimum expected length of the repeat. Default: 25", type=str, default = 25)
parser.add_argument('-M', '--maxlen', help="Maximum expected length of the repeat. Default: 50", type=str, default = 50)

# Extras
parser.add_argument('-v', "--version", action='version', version='%(prog)s '+ versiondisplay)


# Parse resulting arguments
args = parser.parse_args()
# ============================

treksdic = {}
with open(args.treks, 'r') as treks:
	for line in treks: 
		if '>' in line: # Assuming Uniprot database
			db, seqID, fulldescription = line.rstrip("\n").split("|") 
			annotation, taxonomy = fulldescription.split("[") # Extract species name
			taxonomy = taxonomy.rstrip("]")
			
			treksdic[seqID] = (taxonomy, annotation, [])
			newseq = True
		elif 'Length:' in line: #the repeats start
			tabs = line.rstrip("\n").split(" ")
			lenrepeat = int(tabs[1])
			nb = int(tabs[5])
			psim = tabs[13].lstrip('Psim:')
			replen = int(tabs[15].lstrip('Length:'))

			if lenrepeat >= args.minlen and lenrepeat <= args.maxlen:
				treksdic[seqID][2].append([lenrepeat, nb, psim, replen])
				# print(lenrepeat, nb, psim, replen)


for seqID in treksdic.keys():
	taxonomy, annotation, predictions = treksdic[seqID]
	maxregionlen = 0
	maxindex = 0
	if len(predictions) != 0:
		for repindex in range(0, len(predictions)):
			if predictions[repindex][3] > maxregionlen:
				maxregionlen = predictions[repindex][3]
				maxindex = repindex
		
		winner = predictions[maxindex]
		print(f'{seqID}\t{taxonomy}\t{winner[0]}\t{winner[1]}\t{winner[2]}\t{winner[3]}')
	# else: # Report the ones without surviving predicted repeats
	# 	print(f'{seqID}\t{taxonomy}\tNA\tNA\tNA\tNA')



