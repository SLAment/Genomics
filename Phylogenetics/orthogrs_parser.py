#!/usr/bin/env python
# encoding: utf-8

# ================== orthogrs_parser =================
# Script to parse the Orthogroups.csv output file of Orthofinder and manage it
# for the Podospora project.
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/04/08
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
# import sys
import argparse # For the fancy options
import random
# ------------------------------------------------------
version = 1.2
versiondisplay = "{0:.2f}".format(version)

# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Parse Orthogroups.csv for Podospora *", epilog="") # Create the object using class argparse

# Add options
parser.add_argument('orthogrscsv', help="The Orthogroups.csv output file of Orthofinder")
parser.add_argument("--ref", "-r", help="Reference sample. Default: Podan2", default="Podan2")
parser.add_argument("--nugrps", "-n", help="Number of orthologs per species per orthogroup. Default: 1", type=int, default=1)
parser.add_argument("--sample", "-s", help="Sample this number of orthologs. Default: all", type=int, default=float('inf'))
parser.add_argument("--outputdir", "-o", help="Path for output directory", default = ".")
parser.add_argument('--version', "-v", action='version', version='%(prog)s ' + versiondisplay)
parser.add_argument('--verbose', '-b', help="Give some extra information", default=False, action='store_true')

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	orthogrscsvopen = open(args.orthogrscsv, 'r')
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()

# ------------------------------------------------------
# Parse
# ------------------------------------------------------
tabs = [line.rstrip("\n").split("\t") for line in orthogrscsvopen] 			# Read tab file into a list

samples = tabs[0] # Notice the first field is empty # 14

orthogroups = {} # Make a master dictionary 
orthgrlist = [] # keep record in the input order

for line in tabs[1:]:
	# print(line)
	orthogr = line[0]
	orthgrlist.append(orthogr) # Keep record of all the orthogroups

	orthogroups[orthogr] = {} # Make a dictionary for each orthogroup
	for i in range(1,len(samples)): # Ignore the name of the orthogroup
		listgenes = line[i].split(", ")
		if listgenes != ['']: orthogroups[orthogr][samples[i]] = line[i].split(", ") # add nested dictionaries per species (ignore cases without homolog)

# ------------------------------------------------------
# Filter 
# ------------------------------------------------------
niceorthos = []
for ortho in orthgrlist:
	orthogroups[ortho]

	nice = True
	# Does it have all species?
	if list(orthogroups[ortho].keys()) == samples[1:]: 
		for sp in samples[1:]:
			if (len(orthogroups[ortho][sp]) != args.nugrps): nice = False
		if nice: niceorthos.append(ortho) # Append to the final list if nice stayed true

# ------------------------------------------------------
# Sample 
# ------------------------------------------------------
totalnumorthos = len(niceorthos)
if args.sample < float('Inf'):
	samporthos = [] # Sample orthologs
	for x in range(args.sample):
		randnum = random.randint(1, len(niceorthos) + 1) # Generate a random number to use it as index
		samporthos.append(niceorthos[randnum])

	# Replace that full list with the sample:
	niceorthos = samporthos

# ------------------------------------------------------
# Print output
# ------------------------------------------------------
# Print the filtered list of orthogroups
ofile = open(args.outputdir + "/" + "orthogroups_" + str(args.nugrps) + "n.txt", 'w')
for ortho in niceorthos:
	ofile.write(ortho + '\n')

# Recover the name of the reference homolog for each orthogroup
reffile = open(args.outputdir + "/" + args.ref + "_" + str(args.nugrps) + "n.txt", 'w')
for ortho in niceorthos:
	string = ''
	for homolog in orthogroups[ortho][args.ref]:
		string += ''.join(str(homolog)) + "\t"
	reffile.write(string.rstrip("\t") + '\n')

# Print some info
if args.verbose:
	print("Total number of orthogroups: %d" % len(orthgrlist))
	print("Number of orthogroups with %d gene(s) per species: %d" % (args.nugrps, totalnumorthos))
	if args.sample < float('Inf'): print("   Number of sampled orthogroups with %d gene(s) per species: %d" % (args.nugrps, len(niceorthos)))
	print("Output file with orthogroups names: %s" % args.outputdir + "/" + "orthogroups_" + str(args.nugrps) + "n.txt")
	print("Output file with reference names: %s" % args.outputdir + "/" +  args.ref + "_" + str(args.nugrps) + "n.txt")



