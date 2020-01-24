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
import re
# ------------------------------------------------------
version = 1.50
versiondisplay = "{0:.2f}".format(version)

# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Parse Orthogroups.csv for Podospora *", epilog="The result is two files, one with the selected orthogroups, and another with the name of the genes for each orthogroup present in the reference sample (default Podan2).") # Create the object using class argparse

# Add options
parser.add_argument('orthogrscsv', help="The Orthogroups.csv output file of Orthofinder")
parser.add_argument("--sppmap", "-m", help="A species map for the samples in the style of ASTRAL. The filtering is made based on the species, not on the sample (NOT WORKING YET)", type=str)
parser.add_argument("--ref", "-r", help="Reference sample. Default: Podan2", default="Podan2")
parser.add_argument("--nugrps", "-n", help="Number of orthologs per sample per orthogroup. Default: 1", type=int, default=1)
parser.add_argument("--sample", "-s", help="Sample this number of orthologs randomly. Default: all.", type=int, default=float('inf'))
parser.add_argument("--cool", "-c", help="Get all 'cool' orthogroups.", default=False, action='store_true')
parser.add_argument("--coolplus", "-C", help="With cool, print the orthogroups to the screen (for debugging)", default=False, action='store_true')

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

# This patter is Podospora-specific
pagenepattern = re.compile("(P[a-z]+_[\d]_[\d]*)([\.]?[\d]?)_([\w\.]*)") # Exclusive for Podospora annotation, matching genes like Pc_2_8070.2_PcWa133m

# ------------------------------------------------------
# Is there a species map?
# ------------------------------------------------------
if args.sppmap: # Make a dictionary of each species with the corresponding individuals
	sppmp = [line.rstrip("\n").split(":") for line in open(args.sppmap, 'r') ] 			# Read tab file into a list
	sppmapdic = dict([(spp[0], spp[1].split(",")) for spp in sppmp])
	species = list(sppmapdic.keys())

# ------------------------------------------------------
# Parse
# ------------------------------------------------------
tabs = [line.rstrip("\n").split("\t") for line in orthogrscsvopen] 			# Read tab file into a list

samples = tabs[0] # Notice the first field is empty # 14

orthogroups = {} # Make a master dictionary 
orthgrlist = [] # keep record of the input order

for line in tabs[1:]: # Exclude the header
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
orthoswithmissingdata = [] # To keep track of missing data

for ortho in orthgrlist:
	# Filter for the one-to-one orthologs
	nice = True
	complete = True
	samplesinthisortho = list(orthogroups[ortho].keys())

	# --- Does it have all species? ----
	if args.sppmap: 
		if len(samplesinthisortho) < len(species): # It doesn't have all the species for sure
			nice = False
			complete = False
		else: # For the remaining, maybe it does, maybe it doesn't, so loop
			# The orthogroup survives if at least one strain per species has the specified number of orthologs, and none has more than that
			atleastoneguy = [False] * len(species) # At least one individual has to have the requested number of orthologs per species
			for sp in species: 
				intersection = set(sppmapdic[sp]) & set(samplesinthisortho)

				if len(intersection) > 0: # There is an intersection, so there is at least one member of the species
					for sample in intersection: # Check that no sample has more than the requested number of orthologs
						if (len(orthogroups[ortho][sample]) > args.nugrps): # One sample has more than the requested number
							nice = False
						elif len(orthogroups[ortho][sample]) == args.nugrps: # At least one sample has exactly the required number
							atleastoneguy[species.index(sp)] = True		
				else: # If there is no intersection, then at least one of the species is not present in the orthogroup
					nice = False
					complete = False


				# if set(sppmapdic[sp]).issubset(samplesinthisortho): # Ok, so it has all the species with all the individuals, now filter for real
				# 	for sample in sppmapdic[sp]: # Check that no sample has more than the requested number of orthologs
				# 		if (len(orthogroups[ortho][sample]) > args.nugrps): # One sample has more than the requested number
				# 			nice = False
				# 		elif len(orthogroups[ortho][sample]) == args.nugrps:
				# 			atleastoneguy[species.index(sp)] = True
				# 			# print(orthogroups[ortho][sample])
				# # elif len(intersection) > 0: # There is an intersection, so there is at least one member of the species
				# elif len(intersection) == 0: # If there is no intersection, then at least one of the species is not present in the orthogroup
				# 	# print(sppmapdic[sp], samplesinthisortho)
				# 	# print()
				# 	nice = False
				# 	complete = False
				# else: # It's not nice because it's missing individuals, but it still has all species so it's "complete"
				# 	nice = False

			### Useful for debugging
			# if sum(atleastoneguy) == len(species):
			# 	print(atleastoneguy)
			# 	print(species)
			# 	print(orthogroups[ortho])
			# 	print()

		if not complete: orthoswithmissingdata.append(ortho)

	# ---- We treat each sample as a species ----
	else:
		if samplesinthisortho == samples[1:]: 
			for sample in samples[1:]:
				if (len(orthogroups[ortho][sample]) != args.nugrps): nice = False
		else:
			 nice = False
			 orthoswithmissingdata.append(ortho)

	# ---- Finally, check if this ortho passed the test ----
	if args.sppmap:
		if nice and (sum(atleastoneguy) == len(species)):
			niceorthos.append(ortho) 
	elif nice: 
		niceorthos.append(ortho) # Append to the final list if nice stayed true



# ------------------------------------------------------
# Get the cool set of orthologs 
# ------------------------------------------------------
## What makes an ortho cool? 
# - It excludes in-paralogs (duplications within species)
# - It includes groups of paralogs, but not the exact same number per species

if args.cool:	
	if args.verbose: print("Searching for cool paralogs ...")

	coolorthos = []
	singletons = [] # Orthogroups were only one individual has a duplication
	# dupletons = [] # Orthogroups were only one individual has a duplication
	# tripletons = [] # Orthogroups were only one individual has a duplication
	for ortho in orthgrlist:
		cool = True
		samplesinthisortho = list(orthogroups[ortho].keys())

		if not args.sppmap:
			# Does it have all samples?
			if samplesinthisortho == samples[1:]: # yes
				# Make a list of all the genes in this orthogroup
				allgenes = []
				for lista in orthogroups[ortho].values():
					allgenes.extend(lista)

				# Does it have only one gene per species plus only one of the species has two, i.e. singleton inparalogs?
				if len(allgenes) <= (len(samples[1:]) + 1):
					cool = False
					if len(allgenes) == (len(samples[1:]) + 1): # Record the singleton
						singletons.append(ortho)
			else:
				cool = False
		## --- with args.sppmap ---		
		else:
			# At least one individual has to have one ortholog per species
			if ortho in orthoswithmissingdata:
				cool = False
			else: 
				# Is it a singleton?
				atleastoneguy_cool = [False] * len(species) # At least one individual has to have one orthologs per species
				for sp in species:
					intersection = set(sppmapdic[sp]) & set(samplesinthisortho)

					if len(intersection) > 0: # There is an intersection, so there is at least one member of the species
						for sample in intersection: 
							if (len(orthogroups[ortho][sample]) > 1): # It has paralogs
								atleastoneguy_cool[species.index(sp)] = True
					else: # It's missing species
						cool = False

				if sum(atleastoneguy_cool) <= 1: # Then it's a singleton
					cool = False
					if sum(atleastoneguy_cool) == 1: singletons.append(ortho)

				# ## Some extra reporting for me but it's not relevant to the program
				# elif sum(atleastoneguy_cool) == 2:
				# 	# Which species have paralogs
				# 	pairs = [i for (i, v) in zip(species, atleastoneguy_cool) if v]
				# 	if pairs == ['podospora_anserina', 'podospora_comata']:
				# 		print(orthogroups[ortho])

		if cool: coolorthos.append(ortho) # Append to the final list if cool stayed true


	## ---- Deal with the singletons per sample ----
	# Start a dictionary to count how many singletons each sample has
	singledic = dict([(sample, 0) for sample in samples[1:]])
	badannotation_sing = []
	
	for ortho in singletons:
		for sample in list(orthogroups[ortho].keys()): # Who has the extra gene?
			if len(orthogroups[ortho][sample]) > 1:
				# print(ortho, sample, orthogroups[ortho][sample])
				
				# Are they annotation artifacts?
				matchgene1 = pagenepattern.search(orthogroups[ortho][sample][0])
				matchgene2 = pagenepattern.search(orthogroups[ortho][sample][1])

				if matchgene1 and matchgene2: # Are they the same gene but splitted in two?
					if matchgene1.group(1) == matchgene2.group(1): 
						badannotation_sing.append(ortho)
					else: # If they are not, then keep them
						singledic[sample] += 1
						
				else: # I don't know if they are artifacts, so count them
					singledic[sample] += 1
	## ----------------------------------

# print(len(coolorthos))


# ------------------------------------------------------
# Sample 
# ------------------------------------------------------
totalnumorthos = len(niceorthos)
if (args.sample < float('Inf')) and (totalnumorthos > args.sample):
	samporthos = [] # Sample orthologs

	for x in range(args.sample):
		randnum = random.randint(0, totalnumorthos - 1) # Generate a random number to use it as index, endpoints included
		samporthos.append(niceorthos[randnum])

	# Replace that full list with the sample:
	niceorthos = samporthos

# ------------------------------------------------------
# Print output
# ------------------------------------------------------
# Define some printing functions
def writeortholist(ortholist, outputname):
	outfile = open(outputname, 'w')
	for ortho in ortholist:
		outfile.write(ortho + '\n')

def writehomologlist(ortholist, outputname):
	outfile = open(outputname, 'w')
	for ortho in ortholist:
		string = ''
		try: # In case Podan2 (the reference) didn't make it
			for homolog in orthogroups[ortho][args.ref]:
				string += ''.join(str(homolog)) + "\t"
		except:
			string = '?' 
		outfile.write(string.rstrip("\t") + '\n')


# Print the filtered list of orthogroups
writeortholist(niceorthos, args.outputdir + "/" + "orthogroups_" + str(args.nugrps) + "n.txt")
# Recover the name of the reference homolog for each orthogroup
writehomologlist(niceorthos, args.outputdir + "/" + args.ref + "_" + str(args.nugrps) + "n.txt")

# Print the cool orthogroups if pertinent:
if args.cool:
	writeortholist(coolorthos, args.outputdir + "/orthogroups_cool.txt")
	writehomologlist(coolorthos, args.outputdir + "/" + args.ref + "_cool.txt")

	# Also print the genes with the missing data
	writeortholist(orthoswithmissingdata, args.outputdir + "/orthogroups_missing.txt")
	writehomologlist(orthoswithmissingdata, args.outputdir + "/" + args.ref + "_missing.txt")

	if args.coolplus:
		for ortho in coolorthos:
			# Some reporting
			allgenes = []
			for lista in orthogroups[ortho].values():
				allgenes.extend(lista)
			print("** Orthogroup " + ortho + " has " + str(len(allgenes)) + " members **" )


# Print some info
if args.verbose:
	print("Total number of orthogroups: %d" % len(orthgrlist))

	if args.sppmap:
		unit = "species"
		print("Number of orthogroups with %d gene(s) in at least one sample per %s: %d" % (args.nugrps, unit, totalnumorthos))
	else:
		unit = "sample"
		print("Number of orthogroups with %d gene(s) per %s: %d" % (args.nugrps, unit, totalnumorthos))
	
	print("Number of orthogroups with missing genes from at least one %s: %d" % (unit, len(orthoswithmissingdata)) ) # missing data might actually be gene losses 
	if args.sample < float('Inf'): print("   Number of sampled orthogroups with %d gene(s) per sample: %d" % (args.nugrps, len(niceorthos)))
	print("   Output file with orthogroups names: %s" % args.outputdir + "/" + "orthogroups_" + str(args.nugrps) + "n.txt")
	print("   Output file with reference names: %s" % args.outputdir + "/" +  args.ref + "_" + str(args.nugrps) + "n.txt")
	if args.cool: 
		print("Number of orthogroups with no missing data (" + unit + ") that have singletons: " + str(len(singletons)))
		print("   How many singletons are just split genes?", len(badannotation_sing))
		print("   Singletons per sample:", singledic)
		print("Number of cool orthogroups: %d" % (len(coolorthos)))
		print("   Output file with cool orthogroups names: %s" % args.outputdir + "/orthogroups_cool.txt")
		print("   Output file with cool reference names: %s" % args.outputdir + "/" + args.ref + "_cool.txt")

