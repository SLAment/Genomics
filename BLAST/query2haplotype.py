#!/usr/bin/env python3
# encoding: utf-8

# ================== query2haplotype =================
# Script to extract haplotypes of an assembly based on an input query fasta.
# Based on query2hitseq.py and BLAST_extractor.py. It takes as input a fasta
# file, and makes a blastn search to extract the hits sequences from the
# genome. If the --haplo option is used, then it searches for a haplotype
# instead.

# version 2.1 - Removed the use of the Bio.Blast.Applications module (which got depricated) in favor of using subprocess.
# version 2.0 - The script can now do tblastn!
# version 1.92 - added option nocoords to preserve the original contig names if needed
# version 1.9 - added option tophit (but not for --haplo) and fixed the script for python3
# version 1.8 - created option addquery
# version 1.7 - changed the behavior of --self and --noself to check if the output matches the queries
# version 1.6 - added the --makegff option and changed the behaviour of --temp to not include the tailing /
# version 1.5 - added identity argument, and changed the -i flag to -I
# version 1.3 - "from Bio.Alphabet import generic_dna" was removed from BioPython >1.76. See https://biopython.org/wiki/Alphabet
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/05/23, 2021/06/11, 2024/03/14
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import os # For the input name
import sys
from Bio import SeqIO
import subprocess # For the database
import tempfile # To make a temporary file for single sequences
from shutil import rmtree # For removing directories
import argparse # For the fancy options
# ------------------------------------------------------
version = 2.1
versiondisplay = "{0:.2f}".format(version)

# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Extract BLAST hits or haplotypes *", epilog="Warning: the script doesn't work properly if the input sequences have tracks of IUPAC symbols at the beginning and end of the sequences, or '?' symbols anywhere. For python3.") # Create the object using class argparse

# Mandatory options
parser.add_argument('assembly', help="Fasta file to extract from")
parser.add_argument('query', help="Fasta file of the query sequence(s) in nucleotides")
# BLAST options
parser.add_argument("--haplo", "-H", help="Use sequence(s) in query to try to extract haplotypes instead", default=False, action='store_true')
parser.add_argument("--task", "-T", help="Task for the blastn program. Default is 'blastn', other options are 'tblastn', 'blastn-short', 'dc-megablast', 'megablast' or 'vecscreen'", type=str, default="blastn")
parser.add_argument("--evalue", "-e", help="Minimum e-value of BLAST hit to be considered (default 0.00001)", type=float, default=0.00001)
parser.add_argument("--identity", "-i", help="Minimum percentage of identity of BLAST hit to be considered (default 0)", type=float, default=0)
# More filtering options
parser.add_argument("--extrabp", "-f", help="Extra base pairs cut next to the start and end of the BLAST hits or haplotypes (default 0)", type=int, default=0)
parser.add_argument("--minsize", "-s", help="Minimum size of BLAST hit to be considered (default 0 bp)", type=int, default=0)
parser.add_argument("--vicinity", "-c", help="Max distance between 5 and 3 end hits to form a haplotype (default 10000 bp)", type=int, default=10000)
parser.add_argument("--minhaplo", "-m", help="Minimum size of haplotype size (default 0 bp)", type=int, default=0)
# Other useful things
parser.add_argument("--makegff", "-g", help="Print a simple gff3 file with the output hits too", default=False, action='store_true')
parser.add_argument("--addquery", "-q", help="Print the query along with the BLAST results", default=False, action='store_true')
parser.add_argument("--tophit", "-x", help="Slice only the first, top hit of the BLAST search (won't work with --haplo)", default=False, action='store_true')
parser.add_argument("--nocoords", "-0", help="Do not append coordinates of slices into the sequences names", default=False, action='store_true')

# Make a mutualy-exclusive group
selfgroup = parser.add_mutually_exclusive_group()
selfgroup.add_argument("--self", "-N", help="Attempt to report only the hits that are identical to query (won't work if -f is used)", default=False, action='store_true')
selfgroup.add_argument("--noself", "-n", help="Attempt to remove BLAST selfhits (won't work if -f is used)", default=False, action='store_true')

# Extras
parser.add_argument("--blastab", "-b", help="The query is not a fasta, but a BLAST tab file to be used directly instead of doing the whole BLAST", default=False, action='store_true')
parser.add_argument("--seqid", "-I", help="Name of query gene to be extracted from the input multifasta QUERY file. eg. Pa_6_4990", type=str)
parser.add_argument("--threads", "-t", help="Number of threads for BLAST. Default: 1", default="1", type=int) # nargs='+' All, and at least one, argument
parser.add_argument('--temp', '-u', help="Path and directory where temporary files are written in style path/to/dir. Default: working directory", default='.')
parser.add_argument('--clean', '-C', help="Erase the temporary directories and files when finish", default=False, action='store_true')
parser.add_argument('--version', "-v", action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	queryopen = open(args.query, 'r')
	assemblyopen = open(args.assembly, 'r')
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()

# ------------------------------------------------------

# -----------------------------------
## FUNCTIONS
# -----------------------------------
def remove_overlap(ranges):
	""" Simplify a list of ranges; I got it from https://codereview.stackexchange.com/questions/21307/consolidate-list-of-ranges-that-overlap """
	result = []
	current_start = -1
	current_stop = -1 

	for start, stop in sorted(ranges):
		if start > current_stop:
			# this segment starts after the last segment stops
			# just add a new segment
			result.append( (start, stop) )
			current_start, current_stop = start, stop
		else:
			# current_start already guaranteed to be lower
			current_stop = max(current_stop, stop)
			# segments overlap, replace
			result[-1] = (current_start, current_stop) # SLAV: I modified this to update the stop too. (otherwise small ranges contained in a previous larger range will make the current stop smaller than it should be)
	return(result)

# Get a clean name for a file
def namebase(file):
	input_name = os.path.splitext(file)[0] # Remove the extension
	input_name = os.path.basename(input_name) # Remove the path
	return(input_name)


# Make a local database with the given fasta file
def makeBLASTdb(fasta, databasename, dbtype):
	# databasename = name + '_db/' + name + '_db' # Name of the database
	createdb = "makeblastdb -in " + fasta + " -out " + databasename + " -dbtype " + dbtype + " -parse_seqids" # Ugly but compatible with python 2 and 3
	# createdb = f"makeblastdb -in {fasta} -out {databasename} -dbtype {dbtype} -parse_seqids" # BLAST command
	# if args.superverbose: # Print it if necessary
	# 	print("Blast command:", createdb)
	process = subprocess.Popen(createdb.split(), stdout=subprocess.PIPE) # pipe the command to the shell
	stdout, stderr = process.communicate() # run it

# -----------------------------------
# Read input
# -----------------------------------
# Save assembly in memory
records_dict = SeqIO.to_dict(SeqIO.parse(assemblyopen, "fasta"))

# Query
query_dict = SeqIO.to_dict(SeqIO.parse(queryopen, "fasta"))

# If the file comes from RepeatModeler or EDTA, the TE names have a format in the style 
# >TE_00003657_INT#LTR/Copia
# The # breaks the BLAST database.


if args.seqid:
	# Get all sequences that match that string
	matching = [seq for seq in query_dict.keys() if args.seqid in seq]
	if len(matching) != 0: # Make sure there is a match
		# print(matching)
		queryseq = query_dict[matching[0]]
		queryseq.seq.upper()
		# queryseq = query_dict[args.seqid]
	else:
		print("\tNo match for " + args.seqid + "! Exiting ...")
		sys.exit(1)
else:
	queryseq = list(query_dict.keys())

# -----------------------------------
# BLAST query
# -----------------------------------

if not args.blastab:
	# Define the names of the databases
	nameref = namebase(args.assembly)
	nameqry = namebase(args.query)
	databasename = args.temp + "/" + nameref + '_db/' + nameref + '_db'

	# Make the local BLAST databases
	if not os.path.isdir(args.temp + "/" + nameref + '_db/'): # If it exist already, don't bother
		makeBLASTdb(args.assembly, databasename, 'nucl')

	# BLAST
	if args.task in ['blastn-short', 'dc-megablast', 'megablast', 'vecscreen']: 
		args.task = "blastn -task {args.task}"

	if args.seqid: # Only one sequence
		outputhits = args.temp + "/" + queryseq.id + "VS" + nameref + "-" + "hits.tab"

		# Make a temporary file for the fasta file with a single sequence and BLAST it
		with tempfile.NamedTemporaryFile(mode="w", delete=True) as tempquery:
			SeqIO.write(queryseq, tempquery, "fasta")
			tempquery.flush() # important to flush the data because the file might not be closed immediately after writing to it

			blast_command = f"{args.task} -db {databasename} -query {tempquery.name} -out {outputhits} -outfmt 6 -evalue {args.evalue} -num_threads {args.threads} -perc_identity {args.identity}"
			print(blast_command)
			process = subprocess.Popen(blast_command.split(), stdout=subprocess.PIPE) # pipe the command to the shell
			stdout, stderr = process.communicate() # run it

	else: # Multifasta
		outputhits = args.temp + "/" + nameqry + "VS" + nameref + "-" + "hits.tab"

		blast_command = f"{args.task} -db {databasename} -query {args.query} -out {outputhits} -outfmt 6 -evalue {args.evalue} -num_threads {args.threads} -perc_identity {args.identity}"

		process = subprocess.Popen(blast_command.split(), stdout=subprocess.PIPE) # pipe the command to the shell
		stdout, stderr = process.communicate() # run it

	# Read back
	tabs = [line.rstrip("\n").split("\t") for line in open(outputhits, 'r')] 			# Read tab file into a list
else:
	tabs = [line.rstrip("\n").split("\t") for line in open(args.query, 'r')] 			# Read tab file into a list

# -----------------------------------
# Get haplotype
# -----------------------------------

if args.makegff: 
	gfffile = open(args.temp + "/" + nameref + '_vs_' + nameqry + '.gff3', 'w')
	gfffile.write("##gff-version 3\n")
	hitsdic = {}


if args.haplo: # We are looking for entire haplotypes and the blast is only the edges
	chunks = {}
	for hit5 in tabs:
		for hit3 in tabs:
			# Is it the same contig?
			if hit5[1] == hit3[1]:
				hit5start = int(hit5[8])
				hit5end = int(hit5[9])
				hit3start = int(hit3[8])
				hit3end = int(hit3[9])

				maxcoord = max([hit5start, hit5end, hit3start, hit3end])
				mincoord = min([hit5start, hit5end, hit3start, hit3end])
	
				if ((maxcoord - mincoord) <= args.vicinity) and ((maxcoord - mincoord) >= args.minsize):		
					# Add it to dictionary and keep the edges

					if hit5[1] not in chunks.keys():
						chunks[hit5[1]] = [(mincoord, maxcoord)]
					else:
						chunks[hit5[1]].append((mincoord, maxcoord))

	# Slice the hits plus some extra bases on the sides
	slices = []
	chunkcount = 0 # only useful if args.makegff is active
	for ctg in chunks.keys(): # For each contig hit
		hitseq = records_dict[ctg] 	# The actual sequence hit by the query
		for start,end in remove_overlap(chunks[ctg]): # Reduce the list to only the non-overlapping ranges
			chunkcount += 1 # only useful if args.makegff is active
			if (end - start) >= args.minhaplo:
				# Avoid negative numbers
				if start - args.extrabp < 0:
					start_final = 0
				else:
					start_final = start - args.extrabp - 1

				# Stay within the contig
				if end + args.extrabp > len(hitseq.seq):
					end_final = len(hitseq.seq)
				else:
					end_final = end + args.extrabp

				# # Old behaviour
				# if args.noself: # The name of the query and the subject is the same, and the whole of the query sequence
				# 	if (hitseq.id == ctg) and (start == 1) and (end == len(hitseq.seq)): 
				# 		continue # It's a selfhit, so ignore
				# elif args.self:
				# 	if (hitseq.id != ctg) or (start != 1) or (end != len(hitseq.seq)): 
				# 		continue # It's not a selfhit, so ignore

				# Slice it 
				slice = hitseq[start_final:end_final]
				
				if args.nocoords:
					slice.id = hitseq.id 
				else:# Rename it so it's nice
					slice.id = hitseq.id + "_" + str(start_final + 1) + "-" + str(end_final)
					# slice.id = hitseq.id + "_Slice_" + str(start_final + 1) + "-" + str(end_final)
				
				slice.description = ''
				# Save it in the list
				slices.append(slice)

				if args.makegff:
					gfffile.write(f"{hitseq.id}\tBLASTn\tsimilarity\t{start_final + 1}\t{end_final}\t.\t.\t.\tID=haplo{chunkcount};color=#000000;\n")
					# query_dict

else: # The BLAST hits themselves are the haplotypes 
	# Slice the hits plus some extra bases on the sides
	slices = []

	if args.tophit: 
		if len(tabs) > 1: tabs = [tabs[0]] # Assume the first hit is the best one

	for hit in tabs:
		start = min(int(hit[8]), int(hit[9]))
		end = max(int(hit[8]), int(hit[9]))
		hitseq = records_dict[hit[1]] 	# The actual sequence hit by the query

		# Avoid negative numbers
		if start - args.extrabp < 0:
			start_final = 0
		else:
			start_final = start - args.extrabp - 1

		# Stay within the contig
		if end + args.extrabp > len(hitseq.seq):
			end_final = len(hitseq.seq)
		else:
			end_final = end + args.extrabp

		if (abs(start_final - end_final) >= args.minsize):	# (abs(start_final - end_final) <= args.vicinity) and 
			# Slice it 
			slice = hitseq[start_final:end_final]

			if args.nocoords:
				slice.id = hitseq.id 
			else:# Rename it so it's nice
				slice.id = hitseq.id + "_" + str(start_final + 1) + "-" + str(end_final)

			# slice.id = hitseq.id + "_Slice_" + str(start_final + 1) + "-" + str(end_final)
			slice.description = ''
			# Save it in the list
			slices.append(slice)

			if args.makegff:
				# if hitseq.id in hitsdic.keys():
				# 	hitsdic[hitseq.id] += 1
				# else:
				# 	hitsdic[hitseq.id] = 1
				# gfffile.write(f"{hitseq.id}\tBLASTn\tsimilarity\t{start_final + 1}\t{end_final}\t.\t.\t.\tID={hit[0]}_{hitsdic[hitseq.id]};Name={hit[0]};eval={hit[10]};identity={hit[2]};length={len(slice)};color=#000000;\n")
				gfffile.write(f"{hitseq.id}\tBLASTn\tsimilarity\t{start_final + 1}\t{end_final}\t.\t.\t.\tID={hit[0]};Name={hit[0]};eval={hit[10]};identity={hit[2]};length={len(slice)};color=#000000;\n")

if args.noself: # Print only sequences that are different from the queries
	for seq in query_dict:
		for thisslice in slices:
			if thisslice.seq != query_dict[seq].seq: 
				SeqIO.write(thisslice, sys.stdout, "fasta")
elif args.self: # Print only sequences that are 100% identical to the queries
	for seq in query_dict:
		for thisslice in slices:
			if thisslice.seq == query_dict[seq].seq:
				SeqIO.write(thisslice, sys.stdout, "fasta")
else:
	if args.addquery:
		for seq in query_dict.keys():
			SeqIO.write(query_dict[seq], sys.stdout, "fasta")
	# Print the BLAST hits or haplotypes
	SeqIO.write(slices, sys.stdout, "fasta")

# ------------------------------------------------------
# Clean
# ------------------------------------------------------
if args.clean:
	# Remove BLAST databases	
	rmtree(args.temp + "/" + nameref + '_db')
	# Remove the BLAST results
	os.remove(outputhits)
