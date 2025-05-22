#!/usr/bin/env python
# encoding: utf-8

# ================== invertingBEDs =================
# Script to invert the coordinates of a BED file along a scaffold. It works
# with other file types as long as there is a scaffold and at least one
# coordinate column. This is useful when you wish you had done your analysis
# on the reverse-complement sequence of your reference fasta.
# ==================================================
# Sandra Lorena Ament Velasquez
# 2025/05/22
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import argparse  # For the fancy options
import re
import sys  # To exit the script
# ------------------------------------------------------
version = 1.0
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* Script to invert the coordinates of a BED file along a scaffold. *", epilog="The name in the first column of the GFF3 has to be the same as the sequence in the fasta file.")  # Create the object using class argparse

# Add options
parser.add_argument('BED', help="BED or similar file with a scaffold and coordinate columns, tab separated")
parser.add_argument('SCFLENS', help="Support file: a two-columns file with the name of the scaffolds and their total length, tab separated")
parser.add_argument('--fasta', '-f', help="The supporting file is the reference FASTA file used to produce the BED file", default=False, action='store_true')
parser.add_argument("--contig", "-c", help="What column has the contig names? Default: 1", default="1", type=int) # nargs='+' All, and at least one, argument
parser.add_argument("--start", "-s", help="What column has the start coordinates? Default: 2", default="2", type=int) # nargs='+' All, and at least one, argument
parser.add_argument("--end", "-e", help="What column has the end coordinates? Default: 3", default="3", type=int) # nargs='+' All, and at least one, argument

# extras
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	args = parser.parse_args()
	BEDopen = open(args.BED, 'r')
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()

# ------------------------------------------------------

# ---------------------------
### Read the fasta or two-col file
# ---------------------------
records_dict = {}
if args.fasta: # Use the fastas to get the total lenghts of the scaffolds
	from Bio import SeqIO
	for seq_record in SeqIO.parse(args.SCFLENS, "fasta"):
		records_dict[seq_record.id] = len(seq_record)
else:
	for line in open(args.SCFLENS, 'r'):
		scf, length = line.rstrip("\n").split("\t")
		records_dict[scf] = length

# ---------------------------
# Process the BED file
# ---------------------------
newlines = {}
for line in BEDopen:
	if '\t' in line:
		cols = line.rstrip("\n").split("\t")
	else:
		cols = line.rstrip("\n").split()

	scf = cols[args.contig - 1] # The arguments expect a number in base 1, but we're working on base 0 here
	
	try:
		start = int(cols[args.start - 1])
		end = int(cols[args.end - 1])
	except: # This might be a header so there are no numbers in these columns
		sys.stdout.write(line) # print it as is
		continue # move on to the next line

	## Correct the coordinates
	scflen = records_dict[scf]

	if end > start: # depending on the file type, this might be inverted to reflect the orientation of the feature
		invstart = scflen - end + 1 # the +1 is so it starts in base 1 and it's the full window
		invend = scflen - start + 1
	else:
		invstart = scflen - start + 1 # the +1 is so it starts in base 1 and it's the full window
		invend = scflen - end + 1

	if invstart < 0: invstart = 0 # We went over the start, maybe the len of the window is longer than the real contig

	## Reconstruct the line
	newline = ''
	for n in range(0, len(cols)):
		if n == (args.contig-1):
			newline += scf + '\t'
		elif n == (args.start-1):
			newline += str(invstart) + '\t'
		elif n == (args.end-1):
			if args.start != args.end: newline += str(invend) + '\t' # if there is only one column with relevant coordinates, then don't do anything
		else:
			newline += cols[n] + '\t'

	# removing the tailing '\t'
	newline = newline.rstrip('\t') + '\n'

	## Save the line for sorting later
	if scf not in newlines.keys():
		newlines[scf] = [newline]
	else:
		newlines[scf].append(newline)

# ---------------------------
# Print the BED file
# ---------------------------
try: # In case I use `head` with the output # https://stackoverflow.com/questions/26692284/how-to-prevent-brokenpipeerror-when-doing-a-flush-in-python
	for scf in newlines.keys():
		scfranges = newlines[scf]
		scfranges.reverse() # Print in the inverse order:
		for line in scfranges:
			sys.stdout.write(line)

except (BrokenPipeError, IOError): # catch the flush error if head is used
	pass
sys.stderr.close()

