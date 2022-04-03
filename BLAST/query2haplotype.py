#!/usr/bin/env python
# encoding: utf-8

# ================== query2haplotype =================
# Script to extract haplotypes of an assembly based on an input query fasta.
# Based on query2hitseq.py and BLAST_extractor.py. It takes as input a fasta
# file, and makes a blastn search to extract the hits sequences from the
# genome. If the --haplo option is used, then it searches for a haplotype
# instead.

# version 1.5 - added identity argument, and changed the -i flag to -I
# version 1.3 - "from Bio.Alphabet import generic_dna" was removed from BioPython >1.76. See https://biopython.org/wiki/Alphabet
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/05/23, 2021/06/11
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import os # For the input name
import sys
from Bio import SeqIO
from Bio.Blast.Applications import * # NCBI BLASTX wrapper from the Bio.Blast.Applications module
import subprocess # For the database
from shutil import rmtree # For removing directories
import argparse # For the fancy options
# ------------------------------------------------------
version = 1.5
versiondisplay = "{0:.2f}".format(version)

# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Extract BLAST hits or haplotypes *", epilog="Warning: the script doesn't work properly if the input sequences have tracks of IUPAC symbols at the beginning and end of the sequences, or '?' symbols anywhere.") # Create the object using class argparse

# Mandatory options
parser.add_argument('assembly', help="Fasta file to extract from")
parser.add_argument('query', help="Fasta file of the query sequence(s) in nucleotides")
# BLAST options
parser.add_argument("--haplo", "-H", help="Use sequence(s) in query to try to extract haplotypes instead", default=False, action='store_true')
parser.add_argument("--task", "-T", help="Task for the blastn program. Default is 'blastn', other options are 'blastn-short', 'dc-megablast', 'megablast' or 'vecscreen'", type=str, default="blastn")
parser.add_argument("--evalue", "-e", help="Minimum e-value of BLAST hit to be considered (default 0.00001)", type=float, default=0.00001)
parser.add_argument("--identity", "-i", help="Minimum percentage of identity of BLAST hit to be considered (default 0)", type=float, default=0)
# More filtering options
parser.add_argument("--extrabp", "-f", help="Extra base pairs cut next to the start and end of the BLAST hits or haplotypes (default 0)", type=int, default=0)
parser.add_argument("--minsize", "-s", help="Minimum size of BLAST hit to be considered (default 0 bp)", type=int, default=0)
parser.add_argument("--vicinity", "-c", help="Max distance between 5 and 3 end hits to form a haplotype (default 10000 bp)", type=int, default=10000)
parser.add_argument("--minhaplo", "-m", help="Minimum size of haplotype size (default 0 bp)", type=int, default=0)

# Make a mutualy-exclusive group
selfgroup = parser.add_mutually_exclusive_group()
selfgroup.add_argument("--self", "-N", help="Attempt to report only the hits that are identical to query (only with --haplo)", default=False, action='store_true')
selfgroup.add_argument("--noself", "-n", help="Attempt to remove BLAST selfhits (only with --haplo)", default=False, action='store_true')

# Extras
parser.add_argument("--blastab", "-b", help="The query is not a fasta, but a BLAST tab file to be used directly instead of doing the whole BLAST", default=False, action='store_true')
parser.add_argument("--seqid", "-I", help="Name of query gene to be extracted from the input multifasta QUERY file. eg. Pa_6_4990", type=str)
parser.add_argument("--threads", "-t", help="Number of threads for BLAST. Default: 1", default="1", type=int) # nargs='+' All, and at least one, argument
parser.add_argument('--temp', '-u', help="Path and directory where temporary files are written in style path/to/dir/. Default: working directory", default='./')
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
if args.seqid:
	# Get all sequences that match that string
	matching = [seq for seq in query_dict.keys() if args.seqid == seq]
	if len(matching) != 0: # Make sure there is a match
		queryseq = query_dict[args.seqid]
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
	databasename = args.temp + nameref + '_db/' + nameref + '_db'

	# Make the local BLAST databases
	if not os.path.isdir(args.temp + nameref + '_db/'): # If it exist already, don't bother
		makeBLASTdb(args.assembly, databasename, 'nucl')

	# BLAST
	if args.seqid: # Only one sequence
		outputhits = args.temp + queryseq.id + "VS" + nameref + "-" + "hits.tab"
		query_string = '>' + queryseq.id + '\n' + str(queryseq.seq)
		blast_command = NcbiblastnCommandline(cmd='blastn', out=outputhits, outfmt=6, db=databasename, evalue=args.evalue, perc_identity=args.identity, num_threads=args.threads, task=args.task)
		stdout, stderr = blast_command(stdin=query_string)
	else: # Multifasta
		outputhits = args.temp + nameqry + "VS" + nameref + "-" + "hits.tab"
		blast_command = NcbiblastnCommandline(query=args.query, cmd='blastn', out=outputhits, outfmt=6, db=databasename, evalue=args.evalue, perc_identity=args.identity, num_threads=args.threads, task=args.task)
		stdout, stderr = blast_command()

	# Read back
	tabs = [line.rstrip("\n").split("\t") for line in open(outputhits, 'r')] 			# Read tab file into a list
else:
	tabs = [line.rstrip("\n").split("\t") for line in open(args.query, 'r')] 			# Read tab file into a list

# -----------------------------------
# Get haplotype
# -----------------------------------

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

	for ctg in chunks.keys(): # For each contig hit
		hitseq = records_dict[ctg] 	# The actual sequence hit by the query
		cleanchunks = remove_overlap(chunks[ctg])
		for start,end in remove_overlap(chunks[ctg]): # Reduce the list to only the non-overlapping ranges

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

				# if "MarorPlavaka-3" not in hitseq.id: 
				# 	continue

				if args.noself: # The name of the query and the subject is the same, and the whole of the query sequence
					if (hitseq.id == ctg) and (start == 1) and (end == len(hitseq.seq)): 
						continue # It's a selfhit, so ignore
				elif args.self:
					if (hitseq.id != ctg) or (start != 1) or (end != len(hitseq.seq)): 
						continue # It's not a selfhit, so ignore

				# Slice it 
				slice = hitseq[start_final:end_final]
				# Rename it so it's nice
				slice.id = hitseq.id + "_" + str(start_final + 1) + "-" + str(end_final)
				# slice.id = hitseq.id + "_Slice_" + str(start_final + 1) + "-" + str(end_final)
				slice.description = ''
				# Save it in the list
				slices.append(slice)

	SeqIO.write(slices, sys.stdout, "fasta")

else: # The BLAST hits themselves are the haplotypes 
	# Slice the hits plus some extra bases on the sides
	slices = []

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
			# Rename it so it's nice
			slice.id = hitseq.id + "_" + str(start_final + 1) + "-" + str(end_final)
			# slice.id = hitseq.id + "_Slice_" + str(start_final + 1) + "-" + str(end_final)
			slice.description = ''
			# Save it in the list
			slices.append(slice)

	SeqIO.write(slices, sys.stdout, "fasta")

# ------------------------------------------------------
# Clean
# ------------------------------------------------------
if args.clean:
	# Remove BLAST databases	
	rmtree(args.temp + nameref + '_db')
	# Remove the BLAST results
	os.remove(outputhits)
