#!/usr/bin/env python

# ================== DeltaSenseCorrector =================
# Script to detect contigs that are in the wrong direction with respect to a
# reference, using the output (delta file) of the program MUMmer. When there
# are many alignments of one query in one reference, the sense of the longest
# alignment is considered. When the query is mapped to many references, the
# sense of the longest alignment across references is chosen.

# Compared to Version 1: Now it takes options and I can append the names of
# the reference (the chromosomes) to those sequences with unambiguous
# alignment to them (even if it had several alignments within the same
# reference sequence).

# $ /proj/b2015200/Pa_scripts/Alignments/DeltaSenseCorrector.py Dummy.delta dummy.fasta
# $ /proj/b2015200/Pa_scripts/Alignments/DeltaSenseCorrector.py PaWa1m_nc.delta PaWa1m.ctg.fasta
# $ python DeltaSenseCorrector.py /Users/Lorena/Dropbox/PhD_UU/Analyses/Pa_scripts/TestFiles/Pa130p_nc_dummy.delta /Users/Lorena/Dropbox/PhD_UU/Analyses/Pa_scripts/TestFiles/Pa130p.ctg.dummy.fasta -v -a

# version 3 - change of syntax to python3

# ==================================================
# Sandra Lorena Ament Velásquez
# 2015/10/12-20
# +++++++++++++++++++++++++++++++++++++++++++++++++
# version 3
# ------------------------------------------------------
import sys # For reading the input
import os
from Bio import SeqIO # biopython
# from Bio.Alphabet import generic_dna
import argparse # For the fancy options
# ------------------------------------------------------

# ============================
# Make a nice menu for the user
# ============================
version = 3.01
versiondisplay = "{0:.2f}".format(version)

# Make a parser for the options
parser = argparse.ArgumentParser(description="*** DeltaSenseCorrector ***", epilog="Script to detect contigs that are in the wrong direction with respect to a\nreference, using the output (delta file) of the program MUMmer. When there\nare many alignments of one query in one reference, the sense of the longest\nalignment is considered. When the query is mapped to many references, the\nsense of the longest alignment across references is chosen.") # Create the object using class argparse

# Add basic options
parser.add_argument('deltafile', help="Delta file produced with MUMmer")
parser.add_argument('fastafile', help="Fasta file containing the query sequences used with MUMmer")
parser.add_argument('-o', '--outputname', help="Name for the output file including extension (default appends '_fix' to input file name)")

# More
parser.add_argument('-a', '--appendrefname', help="Attach the name of the reference that a query was mapped to (only when there is a single reference) to the query name in the output fasta", default=False, action='store_true')
parser.add_argument('-m', '--onlymapped', help="Print in the new fasta only the sequences that were mapped by MUMmer", default=False, action='store_true')

# Extras
parser.add_argument('-v', "--verbose", default=False, action='store_true')
parser.add_argument('-b', "--version", action='version', version='%(prog)s '+ versiondisplay)


# Parse resulting arguments
args = parser.parse_args()
# ============================


allcontigs = {}
queryrefs = {}
mappedinmanyrefs = []

with open(args.deltafile, 'r') as delta:
    next(delta)									# Header: The files' paths
    next(delta)									# Header: NUCmer or PROmer
    for line in delta:							# The rest of the delta file

        if line[0] == '>':						# if we are dealing with the sequence names
        	refname = line.split(' ')[0][1:]	# get the name of the reference sequence
        	contigname = line.split(' ')[1]		# and the name of the query sequence

        	if contigname in allcontigs.keys(): # If the query sequence is already in our big dictionary
        		allcontigs[contigname].append({refname: [] }) #Append this (expected) unique query sequence name as an internal dictionary
        	else:
	        	allcontigs[contigname] = [{refname: [] }] 	# Else, make an internal dictionary to save all the reference sequences that that query sequence has
        else:
            if ' ' in line: 					# The line with alignment information (with many numbers, not just one)
            	start, end = line.split(' ')[2:4] # Record the coordinates of the query
            	
            	for refdic in allcontigs[contigname]: # For each reference sequence that the query was mapped to
            		if refname in refdic.keys():		# If the current reference sequence is in the contigs
	            		refdic[refname].append((int(start), int(end))) # Save the coordinates of this alignment into a tuple (there can be many alignments per query per reference)

# print(allcontigs)

delta.close()

wrongsense = []
manyrefs = 0
singleref = {} # Record the queries mapped to a single reference

for queryID, references in allcontigs.items():
	# print("**** QueryID:", queryID, "****")
	# print(references)
	start, end = 0, 0

	# Go trough the length of all hits and choose the longest
	currentref = ''
	for refdic in references:
		for refID in refdic:
			for alignment in refdic[refID]:
				lenght = abs(alignment[1] - alignment[0]) 
				oldlenght = abs(end - start)
				if lenght > oldlenght: 					# Is the latest piece longer? else, ignore
					start = alignment[0]				# yes? then these are the new coordinates
					end = alignment[1]
					currentref = refID					# And the new reference
					if start > end:						# Reverse
						senseIsForward = False
					else:								
						senseIsForward = True	
		# print(currentref, start, end, senseIsForward) # currentref is the reference that was accepted as correct
	if senseIsForward == False:					# If the longest piece of this contig is in Reverse sense
		wrongsense.append(queryID)

	if len(references) != 1:	# The query sequences were mapped to many reference sequences
		manyrefs += 1			# Count them
	else:						# The sequences that were mapped to a single reference
		refID = list(references[0].keys())[0] # Get the name of the only reference
		singleref[queryID] = refID 		# Save it in a dictionary

# The name for the output
if args.outputname:
	basenameinput = args.outputname
else:
	basenameinput = os.path.basename(args.fastafile).split(".")[0] + "_fix.fa"

# Give more information if asked
if args.verbose:
	print("** Number of query sequences:", len(allcontigs.keys()), "**")
	print("** Sequences mapped to more than one reference:", manyrefs, "**")
	print("** Sequences reversed complemented:", len(wrongsense), "**")
	if args.appendrefname:
		print("** Number of sequences with the reference name appended:", len(singleref.keys()), "**")

	print("Name of output:", basenameinput)
	if args.onlymapped:
		print("(Only query sequences printed in output fasta file)")

# Open the output file
output_handle = open(basenameinput, "w")

# Write new fasta file correcting the sense of some sequences
for seq_record in SeqIO.parse(args.fastafile, "fasta"):
	if seq_record.id in wrongsense:
		newseq = seq_record.reverse_complement()

		if args.appendrefname and seq_record.id in singleref.keys(): # If asked and the query was mapped to only one reference
			newseq.id = seq_record.id + '_' + singleref[seq_record.id] # append the name of the reference to the query
		else: 
			newseq.id = seq_record.id

		newseq.description = "" # Remove
		# newseq.description = seq_record.description
		SeqIO.write(newseq, output_handle, "fasta")
	else:
		if not args.onlymapped:		# There could be more sequences in the fasta file used that are not present in the delta, so write them too in the output
			SeqIO.write(seq_record, output_handle, "fasta")
		elif seq_record.id in allcontigs.keys():			# Print only the queries that were present in the delta file
			SeqIO.write(seq_record, output_handle, "fasta")

output_handle.close()