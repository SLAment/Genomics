#!/usr/bin/env python
# encoding: utf-8

# ================== totalcovergff =================
# Script to obtained the merged coordinates of all models in gff file (eg.
# from RepeatMasker). Basically it produces a bed file from the gff file with
# overlapping features merged.

# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/06/25
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import argparse  # For the fancy options
import sys  # To exit the script, and to pipe out
import os  # For the input name
# ------------------------------------------------------

version = 2.01
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* Produce a BED file from an input GFF with overlapping features merged *", epilog="The GFF3 has to be sorted beforehand!!!!\nThe name in the first column of the GFF3 has to be the same as the sequence in the fasta file.")  # Create the object using class argparse

# Add options
parser.add_argument('GFF', help="GFF3 file")
parser.add_argument("--fasta", "-f", help="The corresponding assembly to calculate percentage of coverage", type=str)
parser.add_argument("--exclude", "-e", help="Exclude contig(s) in format ctg1,ctg2,ctg3", type=str, default = '')
parser.add_argument("--excludestr", "-E", help="Exclude contig(s) that contain the given string", type=str)

parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	gffopen = open(args.GFF, 'r')
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()
# ------------------------------------------------------

# # Check input file
# # ============================
# try:
#     gff = sys.argv[1]
#     gffopen = open(gff, 'r')
# except:
#     print("Usage: \n$ python " + sys.argv[0] + " sortedannotation.gff > mergedannotation.bed")
#     sys.exit(1)
# # ============================

# ----
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
			result[-1] = (current_start, current_stop) # SLAV: I modified this to update the stop too.

	return(result)
# ----

# ============================
# Read gff 
# ============================
# Save the ranges in a dictionary
gffdic = {} # key: chromosome, value: ([start, end])

for line in gffopen:
	if '#' in line:
		pass
	elif line not in ['\n', '\r\n']: # Ignore empty lines
		cols = line.rstrip("\n").split("\t")

		contig = cols[0]
		start = int(cols[3])
		end = int(cols[4])

		if contig in list(gffdic.keys()):
			gffdic[contig].append([start, end])
		else: # contig is new
			gffdic[contig] = [[start, end]]

# Reduce the overlaps
for ctg in gffdic.keys():
	gffdic[ctg] = remove_overlap(gffdic[ctg])

# ============================
# If provided, read fasta into a dictionary
# ============================
if args.fasta:
	from Bio import SeqIO # For the fasta
	from Bio.Alphabet import generic_dna

	# Read fasta
	fastaopen = open(args.fasta, 'r')
	# This stores in memory
	records_dict = SeqIO.to_dict(SeqIO.parse(fastaopen, "fasta", generic_dna))

	# Exclude contigs if necessary
	if args.excludestr:
		nicectgs = [ctg for ctg in records_dict.keys() if args.excludestr not in ctg]
	else:
		exclulist = (args.exclude).split(",")
		nicectgs = [ctg for ctg in records_dict.keys() if ctg not in exclulist]
	
	# How much is gff covering?
	totalcover = 0
	lenassembly = 0

	sys.stdout.write('Contig\tLen_ctg\tlen_gff\tcoverage\n') # header

	for ctg in nicectgs: # for contigs in the fasta
		thischrtotal = 0
		lenctg = len(records_dict[ctg])
		lenassembly += lenctg

		if ctg in gffdic.keys(): # If there are features for that contig
			for interval in gffdic[ctg]:
				leninterval = abs(interval[1] - interval[0])
				totalcover += leninterval
				thischrtotal += leninterval
			pertg = "{0:.3f}".format((thischrtotal/lenctg)*100)

			sys.stdout.write(f'{ctg}\t{lenctg}\t{thischrtotal}\t{pertg}\n')

	pertg = "{0:.3f}".format((totalcover/lenassembly)*100)
	sys.stdout.write(f'Total\t{lenassembly}\t{totalcover}\t{pertg}\n')

else:
	# Print them
	sys.stdout.write('#Contig\tStart\tEnd\n') # header
	for ctg in gffdic.keys():
		for interval in gffdic[ctg]:
			sys.stdout.write('{}\t{}\t{}\n'.format(ctg, interval[0], interval[1]))
