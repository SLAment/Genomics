#!/usr/bin/env python
# encoding: utf-8

# ================== genometblastn =================
# Script to do a tblastn of a query protein fasta into a genome
# tblastn - Protein Query-Translated Subject BLAST
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020/12/03
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
version = 0
versiondisplay = "{0:.2f}".format(version)

# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Extract BLAST hits or haplotypes *", epilog="") # Create the object using class argparse

# Mandatory options
parser.add_argument('assembly', help="Fasta file to extract from")
parser.add_argument('query', help="Fasta file of the query sequence(s) in nucleotides")
# BLAST options
parser.add_argument("--task", "-T", help="Task for the tblastn program. Options are tblastn (default) or tblastn-fast", type=str, default="tblastn")
parser.add_argument("--evalue", "-e", help="Minimum e-value of BLAST hit to be consider (default 0.00005)", type=float, default=0.00005)
# More filtering options
parser.add_argument("--extrabp", "-f", help="Extra base pairs cut next to the start and end of the BLAST hits (default 0)", type=int, default=0)
parser.add_argument("--minsize", "-s", help="Minimum size of BLAST hit to be consider (default 0 bp)", type=int, default=0)
parser.add_argument("--vicinity", "-c", help="Max distance between 5 and 3 end hits to form a haplotype (default 10000 bp)", type=int, default=10000)
parser.add_argument("--minhaplo", "-m", help="Minimum size of haplotype size (default 0 bp)", type=int, default=0)
# Extras
parser.add_argument("--seqid", "-i", help="Name of query gene to be extracted from the input multifasta QUERY file. eg. Pa_6_4990", type=str)
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
	createdb = f"makeblastdb -in {fasta} -out {databasename} -dbtype {dbtype} -parse_seqids" # BLAST command
	print(createdb)
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

# Define the names of the databases
nameref = namebase(args.assembly)
nameqry = namebase(args.query)
databasename = args.temp + nameref + '_db/' + nameref + '_db'

# Make the local BLAST databases
if not os.path.isdir(args.temp + nameref + '_db/'): # If it exist already, don't bother
	makeBLASTdb(args.assembly, databasename, 'nucl')

print(databasename)
# BLAST
if args.seqid: # Only one sequence
	outputhits = args.temp + queryseq.id + "VS" + nameref + "-" + "hits.tab"
	query_string = '>' + queryseq.id + '\n' + str(queryseq.seq)
	blast_command = NcbitblastnCommandline(cmd='tblastn', out=outputhits, outfmt=6, db=databasename, evalue=args.evalue, num_threads=args.threads, task=args.task)
	stdout, stderr = blast_command(stdin=query_string)
else: # Multifasta
	outputhits = args.temp + nameqry + "VS" + nameref + "-" + "hits.tab"
	blast_command = NcbitblastnCommandline(query=args.query, cmd='tblastn', out=outputhits, outfmt=6, db=databasename, evalue=args.evalue, num_threads=args.threads, task=args.task)
	stdout, stderr = blast_command()

# Read back
tabs = [line.rstrip("\n").split("\t") for line in open(outputhits, 'r')] 			# Read tab file into a list

# -----------------------------------
# Get sequences
# -----------------------------------

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