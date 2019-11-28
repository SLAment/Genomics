#!/usr/bin/env python
# encoding: utf-8

# ================== Stitch BLAST+ hits =================
# Script to stitch together blast hits that should be together, since they are
# greatly overlapping, but for some reason get broken during the BLASTing.

# The output puts a "." for no. of mistmatches and gaps in the stitched hit
# since it's not clear to me how to recover this information.

# The script doesn't handle well cases where there are multiple hits at the
# same position but in opposite directions.
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2016/07/06-07
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import os # For the input name
import sys
# ------------------------------------------------------
version = 3.5
versiondisplay = "{0:.2f}".format(version)

# Input from console
try:
	tabfile = sys.argv[1]
except:
	print("Usage: python " + sys.argv[0] + " tabfile.txt")
	print("Version " + versiondisplay)
	sys.exit(1)

tabopen = open(tabfile, 'r')

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
# ------------------------------------------------------

# ------------------------------------------------------
# Read file into list
# ------------------------------------------------------
unsorted_tabs = [line.rstrip("\n").split("\t") for line in tabopen] 			# Read tab file into a list

# ------------------------------------------------------
# Sort file in case some subjects have multiple hits in different places
# ------------------------------------------------------
# Make a list of all unique queries
queries = []
for tab in unsorted_tabs:
	query = tab[0]
	if query not in queries: queries.extend([query])

# tabs = [] # the future sorted list
THRESHOLD = 2500 # A threshold for the distance between two hits to be distinct loci
stitchedtab = []

# For every query, sort locally
for query in queries:
	subject_dic = {}
	right_subject_order = [] # To keep the right order
	for tab in unsorted_tabs:
		if query == tab[0]: # Just for this query
			subject = tab[1]
			if subject not in subject_dic.keys():
				subject_dic[subject] = [tab]
				right_subject_order.extend([subject]) 
			else:
				subject_dic[subject].append(tab)

	# Sort the resulting list for that particular query, check every subject and sort it locally
	for sub in right_subject_order:
		currentsubject_sorted = sorted(subject_dic[sub], key = lambda x: int(x[9])) # Sort by the subject_start
		# tabs.extend(currentsubject_sorted) # Save it into the final tabs list

		# ------------------------------------------------------
		# Stitch pieces together
		# ------------------------------------------------------
		# print("Subject:", sub)

		maxalign = 0
		queryStarts = []
		queryEnds = []
		subjectStarts = [] 
		subjectEnds = [] 

		for i in range(0, len(currentsubject_sorted)): # Find the main piece
			# print(i, len(currentsubject_sorted), currentsubject_sorted[i] # For debugging)
			# query_id, subject_id, percent_identity, alignment_length, N_mismatches, N_gaps, query_start, query_end, subject_start, subject_end, evalue, bit_score = currentsubject_sorted[i]
			query_start = int(currentsubject_sorted[i][6])
			query_end = int(currentsubject_sorted[i][7])	

			alignment_length = int(currentsubject_sorted[i][3])
			subject_start = int(currentsubject_sorted[i][8])
			subject_end = int(currentsubject_sorted[i][9])

			queryStarts.append(query_start)
			queryEnds.append(query_end)
			subjectStarts.append(subject_start)
			subjectEnds.append(subject_end)

			# Is this last piece one better than the previous ones?
			if alignment_length > maxalign:
				upper_hit = currentsubject_sorted[i]
				maxalign = alignment_length

			# There is only one clean hit
			if len(currentsubject_sorted) == 1: 
				stitchedtab.append(currentsubject_sorted[i])

			# The hit is broken
			elif i < len(currentsubject_sorted) - 1: 
				next_subject_start = int(currentsubject_sorted[i+1][8])

				if abs(subject_end - next_subject_start) > THRESHOLD: # It's probably not part of the same hit # It was subject_start - next_subject_start before
					# So write down the previous one
					# -----------------
					# The new values for the subject
					# -----------------
					# query_id, subject_id, percent_identity, alignment_length, N_mismatches, N_gaps, query_start, query_end, subject_start, subject_end, evalue, bit_score
					new_tab = upper_hit[0:4]

					# These are our new values of the "unbroken" hit, but I'm not sure how to retrieve the no. of gaps and mistmatches
					new_tab.extend(['.', '.', min(queryStarts), max(queryEnds)])

					# subject_start > subject_end for the upper_hit
					if int(upper_hit[8]) > int(upper_hit[9]): # hit is reversed
						new_tab.extend([max(subjectStarts), min(subjectEnds)])
					else:
						new_tab.extend([min(subjectStarts), max(subjectEnds)])
					
					# Let's leave the e-val and Bit score the same as the upper hit
					new_tab.extend(upper_hit[10:])

					stitchedtab.append(new_tab) # Write it in the final output

					# -----------------
					# Reset for the next hit
					# -----------------
					queryStarts = []
					queryEnds = []
					subjectStarts = [] 
					subjectEnds = [] 

					maxalign = int(currentsubject_sorted[i+1][3])
					upper_hit = currentsubject_sorted[i+1]

			else:
				# The last hit in that subject
				# -----------------
				# The new values for the subject
				# -----------------
				new_tab = upper_hit[0:4]

				# These are our new values of the "unbroken" hit, but I'm not sure how to retrieve the no. of gaps and mistmatches
				new_tab.extend(['.', '.', min(queryStarts), max(queryEnds)])

				# subject_start > subject_end for the upper_hit
				if int(upper_hit[8]) > int(upper_hit[9]): # hit is reversed
					new_tab.extend([max(subjectStarts), min(subjectEnds)])
				else:
					new_tab.extend([min(subjectStarts), max(subjectEnds)])
				
				# Let's leave the e-val and Bit score the same as the upper hit
				new_tab.extend(upper_hit[10:])

				stitchedtab.append(new_tab) # Write it in the final output
			
		# print # For separating the subjects

# ------------------------------------------------------
# Print filtered tab file
# ------------------------------------------------------
input_name = os.path.splitext(tabfile)[0]
output_name_tab = input_name+"_stitch.txt"

print_table(stitchedtab, output_name_tab)



