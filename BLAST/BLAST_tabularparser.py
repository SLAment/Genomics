#!/usr/bin/env python
# encoding: utf-8

# ================== Filtering the tabular output file of BLAST+ =================
# Script to parse the tabular output of BLAST+, to produce a tab file that
# satisfy a number of criteria. One of them includes to adjust the parsing to
# the output of Doug's script parseBlastBestHitLoc.pl. The output is printed
# in a new file.

# Based on script BLASTbasicFilter.py done for Lineus project

# TODO: Refine BLAST searches, by merging hits that actually belong to the
# same contig. From Simão et al. 2015 SOM: "Neighbouring highscoring segment
# pairs (HSPs) from the tBLASTn searches are merged if located within 50 Kb of
# each other, thus defining the span of the genomic regions to be evaluated."

# http://pymotw.com/2/argparse/
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2014/07/09-11
# +++++++++++++++++++++++++++++++++++++++++++++++++
# 

# ------------------------------------------------------
import re # For the sorting function
import os # For the input name
import argparse # For the fancy options
from collections import defaultdict #This tells python that the dictionary contains a list so you can freely append things to the value of each key
# ------------------------------------------------------
version = 2.2
versiondisplay = "{0:.3f}".format(version)

# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Parse the output of a BLAST search with different criteria *", epilog="It's VERY important that you verify if the input file has a header or not!!!\nNormally the output of BLAST doesn't have a header, while the output of parseBlastBestHitLoc.pl does.\nIf you're missing the first hit, the header might be the reason!") # Create the object using class argparse

# Add options
parser.add_argument('tabfile', help="TABULAR file(s) obtained after a BLAST search")
parser.add_argument("--outputname", "-o", help="Output name set by the user")
parser.add_argument("--HSPminlen", "-l", help="Minimum HSP length for a hit to be kept (default 0)", default="0", dest="hsp", type=int) # nargs='+' All, and at least one, argument

parser.add_argument("--hitfrac", "-f", help="Use a minimum of 'hit fraction': HSPlength/querylength (default 0)", default="0", type=float) # nargs='+' All, and at least one, argument
parser.add_argument("--identity", "-i", help="Percentage of identity between hit and subject (only for normal BLAST output!)", default="0", type=float) # nargs='+' All, and at least one, argument

parser.add_argument("--evalue", "-e", help="Minimum e-value for keeping a hit (default 10)", type=float, default=10)
parser.add_argument("--Nhits", "-n", help="Number of hits to be recovered after filtering (default all)", type=int)
parser.add_argument("--filterstring", "-y", help="Keep hits with given string (only valid with -s)")
parser.add_argument("--postdoug", "-s", help="The input is the file obtained after using the script parseBlastBestHitLoc.pl", default=False, action='store_true')
parser.add_argument("--headerinput", "-w", help="Turnoff the header for parseBlastBestHitLoc.pl output, but turn it on for normal BLAST", default=True, action='store_false')
parser.add_argument('--verbose', '-b', help="Give some extra information", default=False, action='store_true')
parser.add_argument('--version', "-v", action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	tabopen = open(args.tabfile, 'rU')
except IOError as msg:  # Check that the file exists
    parser.error(str(msg)) 
    parser.print_help()

if args.verbose:
	print 
	print("Input tab file: " + args.tabfile)

# ------------------------------------------------------
def print_table( list_thing, output_name="output.txt", header='' ):
	ofile = open(output_name, 'w')
	if header != '':
		ofile.write((header+'\n'))
	# else:
	# 	print("No header")
	for line in list_thing:
		hitstring = ''
		for tab in line:
			hitstring += ''.join(str(tab)) + "\t"
		ofile.write((hitstring.rstrip("\t") + '\n'))
	ofile.close()

# Got these functions from http://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside by unutbu user
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]
# ------------------------------------------------------

# ------------------------------------------------------
# Read file into list
# ------------------------------------------------------
tabs = [line.rstrip("\n").split("\t") for line in tabopen] 			# Read tab file into a list

badwords = ["hypothetical", "predic", "putative"]



# ------------------------------------------------------
# Does it have a header?
# ------------------------------------------------------
# Depending on the type of intput file 
# Normal Doug (header): T T | -s
if args.postdoug and args.headerinput: 
	headercondition = True
# Normal BLAST F T    | 
elif not args.postdoug and args.headerinput:
	headercondition = False
# BLAST with header: F F    | -w
elif not args.postdoug and not args.headerinput:
	headercondition = True
# Doug without header: T F	| -s -w
elif args.postdoug and not args.headerinput:
	headercondition = False

# The start depends on the type of file
if headercondition:
	startingpoint = 1
	if args.verbose:
		print("** Header active **")
else:
	startingpoint = 0

# ------------------------------------------------------
# Filter the BLAST tab file using a minimum HSP length
# ------------------------------------------------------

nohitcounter = 0
filteredtab = []
if not args.postdoug:
	# heads = "query_id	subject_id	percent_identity	alignment_length	N_mismatches	N_gaps	query_start	query_end	subject_start	subject_end	evalue	bit_score"
	# filteredtab = [line for line in tabs[startingpoint:] if int(line[3]) >= args.hsp] 		# Filter the tab file with hsp_length
	for line in tabs[startingpoint:]:
		if line == ['']: continue # Ignore empty lines
		if int(line[3]) >= args.hsp:
			filteredtab.append(line)
else:
	# heads = "QUERY ACC	QUERY DESCRIPTION	HIT NUM	HIT NAME	HIT DESCRIPTION	SCORE	HIT SIGNIFICANCE	HIT E-VALUE	FRACTION IDENTITY	QUERY LENGTH	HSP LENGTH	QUERY STRAND	QUERY START	QUERY END	HSP STRAND	HSP START	HSP END	NUM GAPS"
	# filteredtab = [line for line in tabs[startingpoint:] if int(line[10]) >= args.hsp] 		# Filter the tab file with hsp_length
	for line in tabs[startingpoint:]:
		if "No hit" in line:
			# print("Query with No hit found")
			nohitcounter += 1
			filteredtab.append(line)
		else:
			if int(line[10]) >= args.hsp:
				filteredtab.append(line)

# Report it
if args.hsp != 0 and args.verbose:
	print("   Minimum HSP length set to %d bp" % args.hsp)

# ------------------------------------------------------
# Make a dictionary with the tab file for each query
# ------------------------------------------------------
hits = defaultdict(list) #This tells python that the dictionary contains a list so you can freely append things to the value of each key
for line in filteredtab:
	hits[line[0]].append(line[1:])

# ------------------------------------------------------
# Sort the names of the queries
# ------------------------------------------------------
sorthitskeys = sorted(hits.keys(), key=natural_keys) 			# Use the function defined in key= to sort

# ------------------------------------------------------
# Keep only some of the hits
# ------------------------------------------------------
hypo = 0

if args.Nhits or args.hitfrac or (args.evalue < 10) or args.filterstring:
	filteredtab = [] 			# Restart the filtered list of hits
	# Establish the position of the e-value
	if not args.postdoug: 		# If BLAST
		eposition = 9
		# HSPlength/querylength
		hspposition = 2
		queryposition = 6
		perident = 1 # Percentage of identity

	else: 							# If post Doug's script
		eposition = 6
		hspposition = 9
		queryposition = 8
		perident = 7 # Percentage of identity
	
	nhits_counter = 0	
	# Do the filtering
	for key in sorthitskeys: # For each query sequence
		# Define a limit of the number of hits recovered
		if args.Nhits:
			limit = args.Nhits
		else:
			limit = len(hits[key])
		
		for hit in hits[key][:]: 	# for each hit of every query 
			if args.postdoug and hit[2] == "No hit": # If you have No hits (and are working with Doug's script)
				newline = hit   
				newline.insert(0, key) 		# Put the name of the query at the beginning
				filteredtab.append(hit) 				# Just append them as they are

			elif float(hit[eposition]) <= args.evalue: # It's below the e-value threshold?
				fraction = float(hit[hspposition])/float(hit[queryposition])
				if fraction >= args.hitfrac: 					# If the hit fraction is above the threshold
					if args.postdoug and args.filterstring: 	# If the we're dealing with Doug's script output and a string to filter was given
						if args.filterstring in hit[3]: 		# Check if the string is present
							newline = hit   			# If yes, take the hit
							newline.insert(0, key) 		# Put the name of the query at the beginning
							filteredtab.append(newline) 	# Add this new line
							nhits_counter += 1
					# else: # If working with normal output or it's Dougs but there was no string to filter with
					elif float(hit[perident]) >= args.identity:
						newline = hit   			# Take the hit
						newline.insert(0, key) 		# Put the name of the query at the beginning
						filteredtab.append(newline) 	# Add this new line
						nhits_counter += 1

				if args.postdoug:
					for word in badwords:
						if word in hit[4]:
							hypo += 1
			# Limit the number of hits recovered
			if nhits_counter >= limit:
				nhits_counter = 0	
				break

# Report some figures if required
if args.verbose:
	if args.Nhits: print("   Max No. of hits kept per query: %d" % args.Nhits)
	if (args.evalue < 10): print("   Minimum e-value set to %f" % args.evalue)
	if args.hitfrac: print("   Minimum hit fraction set to %f" % args.hitfrac)
	if (args.identity > 0): print("   Minimum percentage of identity set to %.2f" % args.identity)

# ------------------------------------------------------
# Print filtered tab file
# ------------------------------------------------------
input_name = os.path.splitext(args.tabfile)[0]

# Get the header out of the first line of the input file
with open(args.tabfile, 'r') as f:
  header = f.readline().rstrip("\n")

# Check if the outputname has been assigned by the user 
if not args.outputname:
	output_name_tab = input_name+"_filtered.txt"
else:
	output_name_tab = args.outputname

# Print the list in a tab format
if headercondition:
	print_table(filteredtab, output_name_tab, header)
else:
	print_table(filteredtab, output_name_tab)

# ------------------------------------------------------
# Some numbers
# ------------------------------------------------------
if args.verbose:
	# Input
	hitsinput = defaultdict(list) #This tells python that the dictionary contains a list so you can freely append things to the value of each key
	for line in tabs[startingpoint:]:
		hitsinput[line[0]].append(line[1:])
	print("\tNumber of input query sequences: " + str(len(hitsinput.keys())))
	print("\tNumber of input hits: " + str(len(tabs) - startingpoint))
	print 

	# Output
	hitsoutput = defaultdict(list) #This tells python that the dictionary contains a list so you can freely append things to the value of each key
	for line in filteredtab:
		hitsoutput[line[0]].append(line[1:])

	print("\tNumber of output query sequences: ", len(hitsoutput.keys()))
	print("\tNumber of output true hits: ", len(filteredtab) - nohitcounter) # True in oposition of having no hit at all but still being in the file
	print("\tNumber of queries without any hits: ", nohitcounter)

	if args.postdoug and (args.Nhits or args.hitfrac or (args.evalue < 10)): 
		print("\tNumber of hits with %s in description: %i" % (badwords, hypo))
	# ------------------------------------------------------
	# The end
	# ------------------------------------------------------
	print("Output tab file: " + output_name_tab)
	print
