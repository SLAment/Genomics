#!/usr/bin/env python
# encoding: utf-8

# ================== collapsegff =================

# Script to collapse the annotation of TEs (from RepeatMasker) into a single
# element if they are too close to each other, and in the same direction.
# This is a super rough approximation, tho, and prone to errors.

# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020/04/22
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import os # For the input name
import sys
import re
# ------------------------------------------------------
version = 1
versiondisplay = "{0:.2f}".format(version)

# Input from console
try:
	gfffile = sys.argv[1]
except:
	print("Usage: python " + sys.argv[0] + " annotation.gff") 
	print("Version " + versiondisplay)
	sys.exit(1)

gffopen = open(gfffile, 'r')
# ------------------------------------------------------

def get_motif(attributes):
	repeatmasker = re.compile('\w*\s"Motif:([\w]+)"\s[\w\s]*')

	matchy = repeatmasker.search(attributes)
	if matchy:
		return matchy.group(1)
	else:
		return ""
# ------------------------------------------------------
# Read file into list
# ------------------------------------------------------
THRESHOLD = 150 # A threshold for the distance between two hits to be distinct loci

# Make a dictionary with the elements of each contig
contigdic = {}
for line in gffopen:
	cols = line.rstrip("\n").split("\t")
	contig = cols[0]

	if contig in contigdic.keys():
		contigdic[contig].append(cols)
	else:
		contigdic[contig] = [cols]

# Now loop through each contig and fuse elements that are too close together
# and in the same orientation
stitchedgff = []

for contig in contigdic.keys():
	sortedlines = sorted(contigdic[contig], key = lambda x: int(x[3])) # Sort by the start of the feature (in case the gff is not sorted)
	# ------------------------------------------------------
	# Stitch pieces together
	# ------------------------------------------------------
	if len(sortedlines) == 1: # There is only one thing annotated in this contig
		stitchedgff.append(sortedlines[0])
	else:
		skipnext = False

		for i in range(0, len(sortedlines)):
			start = int(sortedlines[i][3])
			end = int(sortedlines[i][4])
			sense = sortedlines[i][6]
			motif = get_motif(sortedlines[i][8])


			if skipnext: # This line was already fused with the previous line, so skip it
				skipnext = False
				continue

			if i < len(sortedlines) - 1: # All but the last line in this contig
				nextstart = int(sortedlines[i+1][3])
				nextend = int(sortedlines[i+1][4])
				nextsense = sortedlines[i+1][6]
				nextmotif = get_motif(sortedlines[i+1][8])

				# Fuse them if they are too close, they are the same thing and they are in the same sense
				if (nextstart - end < THRESHOLD) and (sense == nextsense) and (motif == nextmotif): 
					newline = sortedlines[i][0:3]+ [str(start), str(nextend)] + ["."] + sortedlines[i][6:]
					stitchedgff.append(newline)
					skipnext = True

				else: # add it to the final collection
					stitchedgff.append(sortedlines[i])

			elif i == len(sortedlines) - 1: # The last line (will be skipped if it got fused)
				stitchedgff.append(sortedlines[i])

# ------------------------------------------------------
# Print it back into a gff
# ------------------------------------------------------
sys.stdout.write("##gff-version 2\n")
sys.stdout.write(f"#! Gff collapsed from {gfffile} with collapsegff.py v. {versiondisplay}\n")

for line in stitchedgff:
	sys.stdout.write("\t".join(line) + "\n")
