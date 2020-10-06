#!/usr/bin/env python
# encoding: utf-8

# ================== gtfRM2gff =================
# Script to transform the gtf of RepeatMasker into a basic gff3
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020/10/06
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import argparse  # For the fancy options
import datetime
import os # For the input name
import sys # To exit the script 
import re # Regular expressions
# ------------------------------------------------------
version = 1
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* Extract GFF3 features to fasta file *", epilog="The GFF3 has to be sorted beforehand!!!!\nThe name in the first column of the GFF3 has to be the same as the sequence in the fasta file.")  # Create the object using class argparse

# Add options
parser.add_argument('GTF', help="GTF file from RepeatMasker (-gff option in RepeatMasker)")
parser.add_argument("--color", "-c", help="Color in Hex code assigned to the repeats for visualization in IGV (default #cacbd3)", default='#cacbd3')
parser.add_argument("--scolor", "-s", help="Color in Hex code assigned to the simple repeats for visualization in IGV (default #4b5f81)", default='#4b5f81')

# # extras
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	GFFopen = open(args.GTF, 'r')
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()

# ------------------------------------------------------
## Functions

# Functions taken from Joshua Orvis' script convert_aat_btab_to_gff3.py 
next_ids = {'simple_repeat':1} # Dictionary to keep the count
def get_next_id(type, prefix):
    if prefix == None:
        id = type + ".{0:06d}".format( next_ids[type] )
    else:
        id = prefix + "." + type + ".{0:06d}".format( next_ids[type] )
    next_ids[type] += 1
    return(id)

rmpattern = re.compile('Target "Motif:([\w\(\)-]*)"\s([\d]*)\s([\d]*)') # Usually looks like this 'Target "Motif:(TA)n" 1 25'
def attributesmaker(attributes_string):  
	""" Re-shape the gtf attributes to match a gff3 format """
	matchy = rmpattern.search(attributes_string)
	if matchy:
		repeatname = matchy.group(1)	

		# Is this repeat type already in the dictionary (and not a simple repeat)?
		if (")n" not in repeatname) and ("-rich" not in repeatname) and (repeatname not in next_ids.keys()): # If not, add it
			next_ids[repeatname] = 1
		
		# The simple repeats I can put them all together
		if (")n" in repeatname) or ("-rich" in repeatname):
			thisid = get_next_id('simple_repeat', None)
			colour = args.scolor
		else:
			thisid = get_next_id(repeatname, None)
			colour = args.color

		newatts = f'ID={thisid};Name={matchy.group(1)};Note={matchy.group(2)}_{matchy.group(3)};color="{colour}"'
	else:
		print("I can't parse the attributes; are you sure it's a RepeatMasker output?")
		print(attributes_string)
		sys.exit(1)
	return(newatts)
# ------------------------------------------------------


for line in GFFopen: # Without saving in memory
	if '##gff-version' in line: 
		# Add a line to mark the file with this script
		now = datetime.datetime.now()
		metadata = "##GFF3 made out from " + os.path.basename(args.GTF) + " using " + os.path.basename(sys.argv[0]) + " v. " + str(versiondisplay) + ' on ' + str(now) + '\n'
		sys.stdout.write(line)
		sys.stdout.write(metadata)

	elif ('#' in line) or (line == '\n'):
		sys.stdout.write(line)
	else:
		cols = line.rstrip("\n").split("\t")		# break the line into columns defined by the tab

		# Reformat
		newcols = cols
		newcols[2] = "repeat" # Instead of "similarity"
		newcols[8] = attributesmaker(cols[8])
		newline = '\t'.join(newcols) + "\n"

		# Print
		sys.stdout.write(newline)


		# # Print, but also catch the error if I'm piping (see https://stackoverflow.com/questions/26692284/how-to-prevent-brokenpipeerror-when-doing-a-flush-in-python)
		# try:
		# 	sys.stdout.write(newline)
		# except (BrokenPipeError, IOError):
		# 	pass
