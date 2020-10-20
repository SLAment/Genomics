#!/usr/bin/env python
# encoding: utf-8

# ================== GFFSlicer =================
# Script to extract sections of a GFF file while correcting their coordinates,
# making them relative to the start of a given coordinate. The script is meant
# to be used along with fastaSlicer.py, for example.

# The default is to assume there is a single contig in the gff, but if not
# then the contig of insterest can be specified (-c)

# Version 3: change from python2 to python3
# Version 4: change to base 1 and add a few new functions

# TODO: an option to take a specific contig for the gff
# ==================================================
# Sandra Lorena Ament 
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2016/09/14
# +++++++++++++++++++++++++++++++++++++++++++++++++
# $ ./GFFSlicer.py ../TestFiles/test.gff 20000 30349

# ------------------------------------------------------
import sys # To exit the script 
import os # For the input name
import argparse  # For the fancy options
from argparse import RawTextHelpFormatter # To force argparse.ArgumentParser to recognize \n
import datetime
# ------------------------------------------------------
version = 4.02
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* Slice a GFF3 file *", epilog="The GFF3 has to be sorted beforehand!!!!\nThis version can only deal with ONE scaffold at a time.", formatter_class=RawTextHelpFormatter)  # Create the object using class argparse

# Add options
parser.add_argument('GFF', help="GFF3 file sorted beforehand")
parser.add_argument('inputstart', help="The coordinate of the new start, base 1 inclusive!", type=int)
parser.add_argument('inputend', help="The coordinate of the new end, base 1 inclusive! (default inf)", default=float('inf'), nargs='?', type=int) # nargs='?' makes it optional

parser.add_argument("--contig", "-c", help="Name of the contig that should be sliced (otherwise it will assume there is only one contig in the gff)")
parser.add_argument("--outputname", "-o", help="Output name set by the user")
parser.add_argument('--keepcoords', '-k', help="Slice the GFF3 file but without correcting the coordinates", default=False, action='store_true')
parser.add_argument('--invertcoords', '-i', help="Slice the GFF3 file but invert the coordinates (for reverse-complemented sequences); DOES NOT WORK FOR CDS", default=False, action='store_true')

parser.add_argument('--addcoords', '-a', help="Add the value of inputstart to the coordinates instead of slicing (provide any number as end)", default=False, action='store_true')

parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	GFFopen = open(args.GFF, 'r')
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()
# ============================

inputstart = args.inputstart 
inputend = args.inputend
# ---------------------------------
# Name of output
# ---------------------------------
if not args.outputname:
	input_base = os.path.splitext(args.GFF)[0] # Taking out the prefix of the file
	input_name = os.path.basename(input_base) # Remove the path
	newgff = open(input_name + '_Slice_' + str(inputstart) + '_' + str(inputend) +'.gff', "w")
else:
	newgff =  open(args.outputname, "w")
# ---------------------------------

# --------
goodlines = []

for line in GFFopen: # Without saving in memory, but then I can't access the next line so easily
	if '##gff-version 3' in line: 
		goodlines.append(line)
		# Add a line to mark the file with this script
		now = datetime.datetime.now()
		gffslicerstr = "# Sliced from " + os.path.basename(args.GFF) + " using " + os.path.basename(sys.argv[0]) + " v. " + str(versiondisplay) + ' on ' + str(now) + '\n'
		goodlines.append(gffslicerstr)
	elif '#' in line:
		goodlines.append(line)
	elif line == '\n': # If the line is empty
		pass		
	elif '>' in line: # We reached a fasta sequence (assumes the gff stops there)
		break # Stop reading everything.
	else:
		cols = line.rstrip("\n").split("\t")		# break the line into columns defined by the tab

		contig = cols[0]

		if args.contig: # Is this the contig of interest?			
			if contig != args.contig: # If not, ignore the rest
				continue

		currentstart = int(cols[3])
		currentend = int(cols[4])

		# Make a new one to replace values into it
		newcols = cols 

		if args.addcoords: # Append the start to the coordinates instead of slicing
			newstart = currentstart + inputstart
			newend = currentend + inputstart

			newcols[3:5] = str(newstart), str(newend)

			newline = '\t'.join(newcols) + '\n' # Stitch it together as a line
			goodlines.append(newline)

		else:
			# Check that the feature is inside the slice range
			if (currentstart < inputstart):
				continue # If not, ignore this line
			elif (currentend > inputend): # Breaking here will prevent orphan exons that can fit the range of a last gene that doesn't fully
				break
			else: # It it is in range, keep it while fixing the coordinates
				if args.keepcoords: # Leave the coordinates as in the input file
					goodlines.append(line)
				else: # Correct them to match the new range				
					# ----
					if args.invertcoords: # EXPERIMENTAL; Let's assume the gene is sorted
						newstart = (inputend) - currentend 
						newend = (inputend) - currentstart

						# The strand has to be inverted too
						strand = cols[6]
						if strand == "+":
							newcols[6] = "-"

						elif strand == "-":
							newcols[6] = "+"

					else:
					# ----
						newstart = currentstart - inputstart
						newend = currentend - inputstart		

					# Put the new coordinates in 
					newcols[3:5] = str(newstart + 1), str(newend + 1) # The plus one makes sure that it stays in base 1 (inputend is base 0)
					
					newline = '\t'.join(newcols) + '\n' # Stitch it together as a line
					goodlines.append(newline)
		# --------
# Remove the initial left overs of the previous gene
initialbrokengen = []

for line in goodlines:
	if '#' in line:
		continue
	elif 'Parent=' in line:
		initialbrokengen.append(line)
	else:
		break

for bad in initialbrokengen: goodlines.remove(bad)

# --------
# Write the final file
for line in goodlines: newgff.write(line) 

newgff.close()
