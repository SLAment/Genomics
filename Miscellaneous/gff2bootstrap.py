#!/usr/bin/env python
# encoding: utf-8

# ================== gff2bootstrap =================
# Script to get a feeling for the likelihood of a given element to be
# associated with transposable elements (TEs) in a given genome

# Python 3.6+

# Things to keep in mind:
# - I'm assuming there are no Ns in the assembly
# - I'm assuming that almost never I get windows at the edges, so I 
#   don't compensate for smaller windows, so far I just remove them
# - If there are no features for a given contig in the fasta, I assume nothing
# 	was annotated there (absence), rather than thinking the annotation is
# 	missing (missing data)

# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020/04/20
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import argparse  # For the fancy options
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import numpy as np # For random numbers
import time
import sys
# ------------------------------------------------------

version = 1
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* How associated are certain points in the genome with annotated features? *", epilog="The GFF3 has to be sorted beforehand!!!!\nThe name in the first column of the GFF3 has to be the same as the sequence in the fasta file.")  # Create the object using class argparse

# Add options
parser.add_argument('assembly', help="Fasta file with the genomic secuences")
parser.add_argument('focalGFF', help="GFF3 file with focal features to evaluate")
parser.add_argument('GenomeGFF', help="GFF3 file with features of the genome (i.e. other TEs, or genes)")

# Windows
parser.add_argument("--winsize", '-w', help = "Total size of the window (focal point goes in the middle); default 1000 bp", default = 1000, type = int)
parser.add_argument("--minwin", '-m', help = "Minimum fraction of size of edge window to retain it; default 0.75", default = 0.75, type = float)
parser.add_argument("--bootstraps", '-b', help = "Make a distribution of a given number of BOOTSTRAPS; default 100", default = 100, type = int)
parser.add_argument("--noobs", '-N', help = "Print only the bootstraps, not the observed values", default = False, action = 'store_true')
parser.add_argument("--exclude", "-e", help="Exclude contig(s) in format ctg1,ctg2,ctg3", type=str, default = '')

# extras
parser.add_argument('--output', '-o', help="Prefix to append to the output files; eg. 'mysample_'", default = '')
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)
parser.add_argument('--verbose', '-V', help="Print a little bit of extra information", default=False, action='store_true')


try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	fastaopen = open(args.assembly, 'r')
	fGFFopen = open(args.focalGFF, 'r')
	gGFFopen = open(args.GenomeGFF, 'r')
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()

# ------------------------------------------------------

start_time = time.time()

# ----
class Window:
	def __init__(self, coord, size, lenctg): # coord is the center of the window, size is the total length
		self.coord = coord
		self.size = size
		self.finallen = self.size # The final size can change

		if isinstance(coord, list): # If it's a list, order so the order doesn't matter
			lstart = min(coord) # The lists are useful to evaluate the insertion of an element, rather than just a point in the genome
			lend = max(coord)			
		else: 						# If it's just a number (a single point)
			lstart = coord
			lend = coord				

		# Respect the ends of the contig
		if lstart - size/2 < 0: # We're at the edge of the contig
			self.start = 0
			self.finallen = int(self.finallen - (size/2 - lstart)) # Adjust the actual size of the window
		else:
			self.start = int(lstart - size/2)

		if lend + size/2 > lenctg:
			self.end = lenctg
			self.finallen = int(self.finallen - (size/2 - (lenctg - lend))) # Adjust the actual size of the window
		else:
			self.end = int(lend + size/2)

		# Things to add
		self.coverage = None

		# if self.coverage:
		# 	self.percentage = wincoverage/win.finallen

	def __repr__(self):
		if self.coverage:
			return f"Window of size {self.finallen} bp at coordinates {self.coord}; [{self.start}-{self.end}]; coverage {self.coverage}"
		else:
			return f"Window of size {self.finallen} bp at coordinates {self.coord}; [{self.start}-{self.end}]"

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

def get_coverage(gffdic, ctg, win):
	""" How much of that interval is annotated in the gff? """
	wincoverage = 0
	area = 0

	for interval in gffdic[ctg]:
		intervalstart = interval[0]
		intervalend = interval[1]
		if (intervalend < win.start): # the interval is not in the window, but before
			continue
		elif (intervalstart > win.end): # the interval is not in the window, but after (no need to check anymore)
			break
		elif (intervalstart < win.start) and (intervalend > win.start) and (intervalend < win.end): # this interval is partially overlapping, at its end
			area = intervalend - win.start
		elif (intervalstart > win.start) and (intervalend < win.end): # interval fully contained in window
			area = intervalend - intervalstart
		elif (intervalstart > win.start) and (intervalstart < win.end) and (intervalend > win.end): # this interval is partially overlapping, at its beginning
			area = win.end - intervalstart

		wincoverage += area

	return(wincoverage)

# ----

# ---------------------------
### Read the fasta file
# ---------------------------
# # This stores in memory
# records_dict = SeqIO.to_dict(SeqIO.parse(fastaopen, "fasta", generic_dna))

# Exclude contigs if necessary
if args.exclude:
	exclulist = (args.exclude).split(",")
else:
	exclulist = ''

# Make a dictionary with the length of each chromosome	
assemblylens = {seq_record.id:len(seq_record) for seq_record in SeqIO.parse(fastaopen, "fasta") if seq_record.id not in exclulist}
# assemblylens = {seq_record.id:len(seq_record) for seq_record in SeqIO.parse(fastaopen, "fasta")}
	
# ============================
if args.verbose: sys.stdout.write("Reading gff with the assembly features ...\n") 
# ============================
# Save the ranges in a dictionary
gffdic = {} # key: chromosome, value: ([start, end])

for line in gGFFopen:
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

# -----------------------------

# ============================
if args.verbose: sys.stdout.write("Reading gff with the focal features...\n")
# ============================
# Save the ranges in a dictionary
focaldic = {} # key: chromosome, value: ([start, end])

for line in fGFFopen:
	if '#' in line:
		pass
	elif line not in ['\n', '\r\n']: # Ignore empty lines
		cols = line.rstrip("\n").split("\t")

		contig = cols[0]
		start = int(cols[3])
		end = int(cols[4])
		sense = cols[6]

		if contig in list(focaldic.keys()):
			focaldic[contig].append([start, end, sense])
		else: # contig is new
			focaldic[contig] = [[start, end, sense]]


# ============================
## The observed values
# ============================
if not args.noobs: 
	if args.verbose: sys.stdout.write("Writing a file with the observed values...\n")
	obsvalues = open(args.output + 'obsvalues.txt', "w") # Use the current string for the name
	obsvalues.write(f"Contig\tWindow\tStart\tEnd\tDensity\tCoverage\n") # header

obslist = []

### Check the coverage in the window of every focal point
for ctg in focaldic.keys(): # In each chromosome
	for coord in focaldic[ctg]: # Check each window

		lenctg = assemblylens[ctg] # How large is this contig?
		win = Window(coord[:2], args.winsize, lenctg)

		# How much of that interval is annotated in the gff?
		win.coverage = get_coverage(gffdic, ctg, win) # Update window object
		win.percentage = win.coverage/win.finallen

		if win.finallen/args.winsize > args.minwin:
			obslist.append(win.percentage)

		# Print result in a table
		if not args.noobs: obsvalues.write(f"{ctg}\t{win.finallen}\t{win.start}\t{win.end}\t{win.coverage}\t{win.coverage/win.finallen}\n")

	## THE observed value!
	obsmean = sum(obslist)/len(obslist)

sys.stdout.write(f"	... Observed value: {obsmean} (n = {len(obslist)}) at window size of {args.winsize}\n")

# ============================
## The Bootstraps
# ============================

if args.bootstraps > 0:
	sys.stdout.write(f"Calculating a bootstrap distribution ({args.bootstraps} bootstraps)...\n")

	## Choose the chromosome with probability proportional to their sizes
	# This assumes that the values are sorted!
	probs = [prob/sum(assemblylens.values()) for prob in assemblylens.values()]         	# Probabilities proportional to the chromosome sizes

	## Open file for report
	randvalues = open(args.output + f'randvalues_{args.bootstraps}.txt', "w") # Use the current string for the name
	randvalues.write(f"Replicate\tCoverage\n") # header

	# randmean_list = []

	for r in range(args.bootstraps):
		# Some reporting to keep tracking 
		if (r+1)%10 == 0: # +1 because r is in base 0
			sys.stdout.write(f"... random point {r + 1} at window size of {args.winsize}\n")
		
		rnlist = [] # This list will contain the raw coverage points

		countoriginallist = 0 
		while (countoriginallist < len(obslist)): # The while is there to make sure that we always have the same number of random points as the observed data
			rnchr_index = np.random.choice(len(assemblylens.keys()), 1, p=probs, replace=False)[0]	# Choose a chromosome
			rnchr = list(assemblylens.keys())[rnchr_index]											# Get the actual name of that chromosome (contig)

			# Choose the location of the point at random
			rnpoint = np.random.choice(assemblylens[rnchr], 1, replace=False)[0] 

			# Make the window
			lenctg = assemblylens[rnchr] # How large is this contig?
			win = Window([rnpoint,rnpoint], args.winsize, lenctg)

			## How much of that interval is annotated in the gff?
			# Does that contig have any annotation, to begging with?
			if rnchr not in assemblylens.keys(): # no
				win.coverage = 0								# CAREFUL! This assumes absence rather than missing data
			else: # yes
				win.coverage = get_coverage(gffdic, rnchr, win) # Calculate the coverage

			win.percentage = win.coverage/win.finallen

			if win.finallen/args.winsize > args.minwin: # This window is accepted
				rnlist.append(win.percentage)
				countoriginallist += 1

		randmean = sum(rnlist)/len(rnlist)
		
		# Print the bootstrap into a file
		randvalues.write(f'{r + 1}\t{randmean}\n')
		# randmean_list.append(randmean)

	randvalues.close()

# ## Rudimentary plot
# from matplotlib import pyplot as plt

# # Compute frequency and bins
# plt.hist(randmean_list, bins=50)
# plt.gca().set(title='Frequency Histogram', ylabel='Frequency')
# plt.show()

if args.noobs and (args.bootstraps < 1): print("Nothing done because -N and 0 bootstraps :P")

sys.stdout.write("Done!\n")
sys.stdout.write("--- {0:.2f} seconds ---\n".format(time.time() - start_time))
