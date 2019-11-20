#!/usr/bin/env python
# encoding: utf-8

# ================== GFFSlicer =================
# Script to extract sections of a GFF file while correcting their coordinates,
# making them relative to the start of a given coordinate. The script is meant
# to be used along with fastaSlicer.py, for example.

# I didn't consider the issue with the contig in turn. As it is, GFFSlicer.py
# can only work with GFF files of a single chromosome. Further improvement
# could be to make an option to extend the coordinates of a GFF file, instead
# of reducing them.

# Version 3: change from python2 to python3

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
version = 4.00
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* Slice a GFF3 file *", epilog="The GFF3 has to be sorted beforehand!!!!\nThis version can only deal with ONE scaffold per GFF file.", formatter_class=RawTextHelpFormatter)  # Create the object using class argparse

# Add options
parser.add_argument('GFF', help="GFF3 file sorted beforehand")
parser.add_argument('inputstart', help="The coordinate of the new start, base 1 inclusive!", type=int)
parser.add_argument('inputend', help="The coordinate of the new end, base 1 inclusive! (default inf)", default=float('inf'), nargs='?', type=int) # nargs='?' makes it optional

parser.add_argument("--outputname", "-o", help="Output name set by the user")
parser.add_argument('--keepcoords', '-k', help="Slice the GFF3 file but without correcting the coordinates", default=False, action='store_true')

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

## Make the coordinates base 0 (input should be base 1)
if args.inputstart > args.inputend: # Inverted coordinates
	inputstart = args.inputend - 1
	inputend = args.inputstart - 1
else:
	inputstart = args.inputstart - 1
	inputend = args.inputend - 1

# ---------------------------------
# Name of output
# ---------------------------------
if not args.outputname:
	input_base = os.path.splitext(args.GFF)[0] # Taking out the prefix of the file
	input_name = os.path.basename(input_base) # Remove the path
	newgff = open(input_name + '_Slice_' + str(inputstart + 1) + '_' + str(inputend + 1) +'.gff', "w")
else:
	newgff =  open(args.outputname, "w")
# ---------------------------------

# --------
goodlines = []

for line in GFFopen:
	# print(line)
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

		currentstart = int(cols[3])
		currentend = int(cols[4])

		if args.addcoords: # Append the start to the coordinates instead of slicing
			newstart = currentstart + inputstart
			newend = currentend + inputstart

			newcols = cols # Put the new coordinates in
			newcols[3:5] = str(newstart), str(newend)

			newline = '\t'.join(newcols) + '\n' # Stitch it together as a line
			goodlines.append(newline)

		else:
			# Check that the feature is inside the slice range
			if (currentstart <= inputstart): # inputstart is Base 0
				continue # If not, ignore this line
			elif (currentend >= inputend): # Breaking here will prevent orphan exons that can fit the range of a last gene that in totally doesn't
				break
			else: # It it is in range, keep it while fixing the coordinates
				if args.keepcoords: # Leave the coordinates as in the input file
					goodlines.append(line)
				else: # Correct them to match the new range
					newstart = currentstart - inputstart
					newend = currentend - inputstart

					newcols = cols # Put the new coordinates in
					newcols[3:5] = str(newstart), str(newend)

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
