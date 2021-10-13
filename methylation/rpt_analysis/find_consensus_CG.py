### Background info
	### bismarkMultimappingLocations has an option called --limit x that reports up to x multimapping locations
	### for the top strand and up to x location for the bottom (up to 2x location in total)
	### The top/bot strand is the reverse complement of the other

### Algorithm
	### Get *ambig.bam locations of a read
	### Get reference sequence of these locations (bedtools getfasta)
		### Make reverse complement top or bot strand locations so that all sequences
		### face the same direction
	### Run MUSCLE on these reference sequences
	### Check if there are CGs that are conserved in all multimapped locations
	### Calls methylation on CGs

### Requirements
	### NEEDS *AMBIG.BAM & REFERENCE FASTA
		### READS W/ SAME NAME NEED TO BE GROUPED/SORTED
	### NEEDS FASTQ FILE THAT *AMBIG.BAM WAS GENERATED FROM
	### This script only takes bismarkMultimappingLocations *ambig.bam output
	### requires pybedtools
	### requires commandline muscle
	### requires pybam



import pysam
import pybedtools
import argparse
import re
import os
import gzip
from collections import defaultdict
from sys import argv
import csv
from string import maketrans
from subprocess import call
from difflib import SequenceMatcher
import string

################################################################################################################################################***

# funcation to find consensus CpGs in MUSCLE output (clw format)
# Outputs the relative location of every consensus CpG on the string
# A consensus CpG is one that is on the same relative location of the reference for every location a read tried to multimap to
def findconsensus(filename):
	# Lines are separated by \n. Make lines into 1 combined line
	strings = defaultdict(lambda: "")
	with open(filename, 'r') as f:
		for n,line in enumerate(f):
			if n > 0:
				line = line.strip().split()
				if line == []:
					continue
				elif "*" in ''.join(line):
					continue
				else:
					key = line[0].split('_')[0]
					strings[key] = strings[key] + line[1]
	f.close()
	cgloc = {}

	# What are the relative locations of CGs that are shared between reference locations?
	for i in strings.keys():
		cgloc[i] = set([m.start() for m in re.finditer('CG', strings[i])])

	# This list will contain the relative pos of shared CpGs
	consensus = reduce(lambda s1, s2: s1 & s2, cgloc.values())

	# The issue with this is that MUSCLE will insert dashes that shift the actual reference 
	# position of the CpG if not accounted for
	# In addition, some reads are the reverse complement, so the ref position of the CpG is "inverted"
	# Declare new dictionary with "adjusted" relative consensus positions
	adjustpos = defaultdict(lambda: ())
	if consensus != 0:
		consensus = sorted(consensus)
		for keys in strings.keys():
			seq = strings[keys]
			if keys.endswith('+)'): 
				keys = keys.split(':')
				for cgpos in consensus:
					# I only care about how many '-'s are to the left of a CpG as deleting only them affects position
					lenleft = len([m.start() for m in re.finditer('-', seq[:cgpos])])
					adjustpos[':'.join(keys[0:3])] = adjustpos[':'.join(keys[0:3])] + (cgpos - lenleft,)
			# Deal with rverse comp
			elif keys.endswith('-)'): 
				keys = keys.split(':') 
				for cgpos in consensus:
					length = int(keys[2]) - int(keys[1])
					lenleft = len([m.start() for m in re.finditer('-', seq[:cgpos])])
					tmppos = (cgpos-lenleft)
					# This flips the rev comp position
					adjustpos[':'.join(keys[0:3])] = adjustpos[':'.join(keys[0:3])] + (length - 2- tmppos,)
			
	# Return a tuple that contains consensus and adjustpos
	return (consensus, adjustpos,)

################################################################################################################################################***

# Function to deal with insertions and deletions in reads. A key part of this code is translating the relative position of CpGs on a reference sequence to the position on the
# bisulfite converted read (see findconsensus above). However the read can have indels which change the expected position(s). This function returns a value that reflects 
# how much these positions have changed due to indels. Needs cigar tuple and the reference CpG relative postion (index)
def indel(cigartuple, index):
	# Initally, shifted amount (itr) is 0. flagpos is the cigar position(s) in the cigar string. It starts at 1 so we need
	# to substract 1
	itr = 0
	flagpos = -1
	for flag in cigartuple:
		flagpos += int(flag[1])
		if index >= flagpos and int(flag[0]) != 0:
			if int(flag[0]) == 2: 
			#deletion
				itr -= int(flag[1])
				flagpos += int(flag[1])

			elif flag[0] == 1: 
				#insertion
				itr += int(flag[1])
				flagpos -= int(flag[1])
		else:	
			continue
	return itr

################################################################################################################################################***

# Call methylation. 'CG' is methylated, 'TG' or 'CA' is unmethylated
# Takes a sequence, index and itr (see above)
def callmeth(seq, index, itr):
	if seq[index+itr:index+itr+2] == 'CG':
		meth = 1
	elif seq[index+itr:index+itr+2] == 'TG' or seq[index+itr:index+itr+2] == 'CA':
		meth = -1
	else:
		if len(seq[index+itr:index+itr+2]) > 0:
			if len(seq[index+itr:index+itr+2]) < 2:
				meth = 'CUT_OFF'
			else:
				meth = 'MUTATION_' + seq[index+itr:index+itr+2]
		
		else:
			meth = "NOT_COV"

	return meth

################################################################################################################################################***

# Reverse complement
def revcomp(seq):
	trans = string.maketrans('ATGC','TACG')
	return seq.translate(trans)[::-1]

################################################################################################################################################***

def main(args):

	# Parse arguments
	bam = args.bam
	limit = args.limit
	fasta = args.fasta
	poption = args.poption
	fastq = args.fastq

	#############################################################################

	# BismarkMultimappingPositions output .ambig files do NOT agree with fastq sequence. Methylation calls
	# have been erased. Need to obtain the fastq sequence and replace the SAM sequence
	# To do this, I create another dictionary that holds sequence info for every read name. 
	# Make sure these fastq names contain the string '_1:' or '_2:' that denotes which mate the read is

	fastq_seq = {}
	i = 1


	with gzip.open(fastq, 'rb') as f:
		for line in f:
			line = line.strip()
			if i == 1:
				name = line.lstrip('@')
			elif i == 2:
				fastq_seq[name] = line
			elif i ==4:
				i = 0
			i += 1
	print "Sucessfully loaded fastq. \n"

	#############################################################################
	# Parse thru bam
	# I ASSUME that every unique read has at least 2 multimapping locations 
	# Bisamrk ambig reports alignments to the top and bottom strand. However, it doesnt specify
	# which stand the sequence aligned to.
	# Therefore, the first sequence the script sees will be the "top" (+) and all sequences for a read name that is the reverse comp will be aligned to "bottom" (-)

	# Create dictionary to hold read name and positional information as well as sequence orientation info (+-)
	# structure is {read name: set(positions)}
	orig_reads = defaultdict(lambda: ())

	# Create dictionary that holds sequence and cigar info. The reason I am not including this info
	# In the above dics is because that info will go into a bed file and Id rather not split the string
	seq_orig_reads = {}

	with pysam.AlignmentFile(bam, "rb") as f:
		for line_num,line in enumerate(f):
			if line_num == 0:
				orig_seq = line.query_sequence
				tmp_name = line.qname
			current_name = line.qname
			#print current_name
			start = str(line.pos)
			stop = str(int(start) + int(line.query_alignment_length) + 1)
			chrm = line.reference_name
			current_seq = line.query_sequence
			# TEST
			if len(fastq_seq[current_name]) != len(current_seq):
				print "BAM and FASTQ sequence i not the same length"
				quit()

			# choose which dic to add to. Always add to orig_reads 1st due to assumption
			# 
			# To make things even more complex, the reads in BAM and FASTQ arent nessessarily in the same orientation
			# I add a tag to the end of the keys of the strings dictionary that denotes whether or not 
			# I check this by comparing the first 15 char and last 15 revcomp char in the fastq string with the bam sequence
				# If there is more similarity b/t last last 15 revcomp char in the fastq string with the bam sequence
				# We take the reverse complement of the fastq sequence
		
			if current_name == tmp_name:
				fwdratio = SequenceMatcher(None, fastq_seq[current_name][0:15], current_seq[0:15]).ratio()
				rvratio = SequenceMatcher(None, revcomp(fastq_seq[current_name])[0:15], current_seq[0:15]).ratio()
				if current_seq == orig_seq:
					if fwdratio > rvratio:
						seq_orig_reads["%s:%s:%s:%s" % (current_name, chrm, start, stop)] = (fastq_seq[current_name], line.cigartuples,)
					elif fwdratio < rvratio:
						seq_orig_reads["%s:%s:%s:%s" % (current_name, chrm, start, stop)] = (revcomp(fastq_seq[current_name]), line.cigartuples,)
					else: # If ratios equal, try again with larger string
						fwdratio = SequenceMatcher(None, fastq_seq[current_name], current_seq).ratio()
						rvratio = SequenceMatcher(None, revcomp(fastq_seq[current_name]), current_seq).ratio()
						if fwdratio >= rvratio:
							seq_orig_reads["%s:%s:%s:%s" % (current_name, chrm, start, stop)] = (fastq_seq[current_name], line.cigartuples,)
						elif fwdratio < rvratio:
							seq_orig_reads["%s:%s:%s:%s" % (current_name, chrm, start, stop)] = (revcomp(fastq_seq[current_name]), line.cigartuples,)

					tmplist = [chrm, start, stop, '%s:%s:%s:(+)_%s' % (chrm, start, stop, current_name), '1', '+']
					orig_reads[current_name] = orig_reads[current_name] + (tmplist,)
				else:
					if fwdratio > rvratio:
						seq_orig_reads["%s:%s:%s:%s" % (current_name, chrm, start, stop)] = (fastq_seq[current_name], line.cigartuples,)
					elif fwdratio < rvratio:
						seq_orig_reads["%s:%s:%s:%s" % (current_name, chrm, start, stop)] = (revcomp(fastq_seq[current_name]), line.cigartuples,)
					else:
						fwdratio = SequenceMatcher(None, fastq_seq[current_name], current_seq).ratio()
						rvratio = SequenceMatcher(None, revcomp(fastq_seq[current_name]), current_seq).ratio()
						if fwdratio >= rvratio:
							seq_orig_reads["%s:%s:%s:%s" % (current_name, chrm, start, stop)] = (fastq_seq[current_name], line.cigartuples,)
						elif fwdratio < rvratio:
							seq_orig_reads["%s:%s:%s:%s" % (current_name, chrm, start, stop)] = (revcomp(fastq_seq[current_name]), line.cigartuples,)

					tmplist = [chrm, start, stop, '%s:%s:%s:(-)_%s' % (chrm, start, stop, current_name), '1', '-']
					orig_reads[current_name] = orig_reads[current_name] + (tmplist,)
			else:
				orig_seq = current_seq
				tmp_name = current_name
				fwdratio = SequenceMatcher(None, fastq_seq[current_name][0:15], current_seq[0:15]).ratio()
				rvratio = SequenceMatcher(None, revcomp(fastq_seq[current_name])[0:15], current_seq[0:15]).ratio()
				if fwdratio > rvratio:
					seq_orig_reads["%s:%s:%s:%s" % (current_name, chrm, start, stop)] = (fastq_seq[current_name], line.cigartuples,)
				elif fwdratio < rvratio:
					seq_orig_reads["%s:%s:%s:%s" % (current_name, chrm, start, stop)] = (revcomp(fastq_seq[current_name]), line.cigartuples,)
	
				else:
					fwdratio = SequenceMatcher(None, fastq_seq[current_name][0:15], current_seq[0:15]).ratio()
					rvratio = SequenceMatcher(None, revcomp(fastq_seq[current_name])[0:15], current_seq[0:15]).ratio()
					if fwdratio >= rvratio:
						seq_orig_reads["%s:%s:%s:%s" % (current_name, chrm, start, stop)] = (fastq_seq[current_name], line.cigartuples,)
					elif fwdratio < rvratio:
						seq_orig_reads["%s:%s:%s:%s" % (current_name, chrm, start, stop)] = (revcomp(fastq_seq[current_name]), line.cigartuples,)

				tmplist = [chrm, start, stop, '%s:%s:%s:(+)_%s' % (chrm, start, stop, current_name), '1', '+']
				orig_reads[current_name] = orig_reads[current_name] + (tmplist,)
	f.close()
	del fastq_seq
	
	print "Successfully loaded %s. Corresponding sequences in fastq exist." % (bam)

	#############################################################################
	# Delete reads with more multimapping locations than the limit
	if limit > 1:
		for name in orig_reads.keys():
			if (len(orig_reads[name])) > limit:
				orig_reads.pop(name)

	#############################################################################
	# only output consensus CpGs from reads whose mate is multimapped near it
	# parse thru dictionaries and only return sequences with mates
	
	if poption != 'no':
		# Create 2 sets for mates
		r1set = set()
		r2set = set()
		
		# Place read names into sets
		for name in orig_reads.keys():
			if '_1:N:' in name:
				r1set.add(name)
			elif '_2:N:' in name:
				r2set.add(name)

		# Turn r2 names into r1 names and intersect sets to see which reads have mates
		r2tmp = set(map(lambda x: x.replace('_2:N:', '_1:N:', 1), list(r2set)))
		r1mates = r1set & r2tmp
		del (r1set, r2set, r2tmp)
		r2mates = set(map(lambda x: x.replace('_1:N:', '_2:N:', 1), list(r1mates)))
		
		# Delete all reads without mate
		nomate = set(orig_reads.keys()) - r1mates
		nomate = nomate - r2mates
		for name in nomate:
			orig_reads.pop(name)

		# Now only return locations where r1 is multimapped near r2
		# Declare new dictionary to hold only reads with mates near it
		new_reads = defaultdict(lambda: ())

		for r1 in r1mates:
			r2 = r1.replace('_1:N:', '_2:N:', 1)
			for line1 in orig_reads[r1]:
				for line2 in orig_reads[r2]:
					# Declare chrm of r1 and r2
					c1 = line1[0]
					c2 = line2[0]
					if c1 == c2:
						start1 = int(line1[1])
						start2 = int(line2[1])
						if (start2-305) <= start1 and start1 < (start2+305):
							new_reads[r1] = new_reads[r1] + (line1,)
							new_reads[r2] = new_reads[r2] + (line2,)

					else:
						continue
		# Replace orig_reads with new_reads
		del orig_reads
		orig_reads = new_reads
		
		print "Successfully filtered out reads without mates in vicinity. \n"

	#############################################################################	

	# Turn dictionaries into bedfiles into reference fasta
	# MUSCLE fasta to see if there are shared CGs
	# String to hold bed file info
	tmpfasta = bam + ".TMP.fasta.pyes"
	tmpmuscle = bam + ".TMP.musc.pyes"
	tmpoutput = bam + ".ConsensusCpG.bed.tmp.pyes"
	output = bam + ".ConsensusCpG.pyes.bed"

	# Loop thru each unique read name and turn into tmporary fasta, and run muscle many times
	# Open output file
	print "Done processing. Now running MUSCLE\n"

	lengthdic = len(orig_reads.keys())
	with open(tmpoutput, 'w') as o:
		o.write("#chr\tstart\tstop\treadLocation_Name\tmethylation\n")
		for progress, name in enumerate(orig_reads.keys()):
			bedstring = """"""
			count = 0
			# Print bed locations to string (bedstring)
			for bed in orig_reads[name]:
				if bed[0].startswith('chr'):
					bedstring = bedstring + "%s\n" % ('\t'.join(bed))
					count += 1
			# continues with next iteration if not enough multimapping locations
			if count < 2:
				continue

			# Run bedtools getfasta (.sequence). The output is now a FASTA file
			pybedobject = pybedtools.BedTool(bedstring, from_string = True)
			pybedobject = pybedobject.sequence(fi=fasta, s=True, name=True)

			# attempt to make code faster by check 1st sequence for cg
			# if no cg, continue with next iteration of loop
			if 'CG' not in open(pybedobject.seqfn).read()[:500].splitlines()[1].upper():
				continue

			with open(tmpfasta, 'w') as f:
				f.write(open(pybedobject.seqfn).read())
			f.close()

			# Now run MUSCLE on FASTA file
			tmp = call(['muscle', '-in', tmpfasta, '-out', tmpmuscle, '-maxiters', '1', '-diags1', '-clw', '-quiet'] )
		
			# Find consensus CpGs and output locations as BED file
			returntuple = findconsensus(tmpmuscle)
			consensusloc = returntuple[0]
			adjustposloc = returntuple[1]

			if len(consensusloc) == 0:
				continue

			# Output as bed file and also output methylation info

			for bed in orig_reads[name]:
				if bed[0].startswith('chr'):
					#bed = bed.strip().split('\t')
					bedpos = ':'.join(bed[0:3])
					for keys in adjustposloc.keys():
						if bedpos == keys:
							for index in adjustposloc[keys]:
								index = int(index)
								start = int(bed[1]) + index
								stop = start + 1
								# Check methylation
								seq_index = name + ':' + bedpos
								# Deal with insertions and deletions
								cigartuple = seq_orig_reads[seq_index][1]
								if len(cigartuple) > 1:
									itr = indel(cigartuple, index)
								else:
									itr = 0

								# Call methylation
								meth = callmeth(seq_orig_reads[seq_index][0], index, itr)

								o.write("%s\t%s\t%s\t%s\t%s\n" % (bed[0], start, stop, bed[3], meth))
									
											
							break
				else:
					pass
			# Output progress
			if progress == int(lengthdic) / 4:
				print ".. 25%\n"
			if progress == int(lengthdic) / 2:
				print ".. 50%\n"
			if progress == (int(lengthdic) *3) / 4:
				print ".. 75%\n"
		o.close

	
	# Remove tmp files
	os.remove(tmpfasta)
	os.remove(tmpmuscle)

	print "100%\nDone finding consensus CpGs and calling methylation.\n"

	#############################################################################

	# Now characterize methylation (calculate coverage and percent methylation for unique CpG locations)
	# cpgdic[positions] = (coverage, numCs)
	cpgdic = defaultdict(lambda: [0,0])

	outfile = open(output, 'w')
	outfile.write("#chr\tstart\tcov\tnumCs\t5mCratio\n")
	with open(tmpoutput, 'r') as f:
		for line in f:
			line = line.split('\t')
			try:
				loc = "%s\t%s\t%s" % (line[0],int(line[1],int(line[2]))
				cpgdic[loc][0] += 1
				if int(line[4]) == 1:
					cpgdic[loc][1] += 1
			except:
				pass
		for l in sorted(cpgdic.keys()):
			outfile.write("%s\t%s\t%s\n" % (l, '\t'.join(map(str, cpgdic[l])), (float(cpgdic[l][1]) / float(cpgdic[l][0]))))


		f.close()

	print "Done characterzing methylation on each unique CpG \n\nWorkflow complete!"

################################################################################################################################################***

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='See comments in code')
        parser.add_argument('-b', dest='bam', type=str, help='*ambig_bam file from bismarkMultimappingLocations')
        parser.add_argument('-f', dest='fasta', type=str, help='Reference fasta name')
	parser.add_argument('-q', dest='fastq', type=str, help='Reference fastq name')
	parser.add_argument('-l', dest='limit', default = 0, type=int, help='Optional: Limit on the # of multimapped locations a read can have (takes integer and run time is longer). If a read has more locations than the specified limit, it is NOT analyzed')
	parser.add_argument('-p', dest='poption', default = 'no', type=str, help='Optional: type "yes" if paired-end data and you would like to only output consensus CpGs from reads whose mate is multimapped near it')
	args=parser.parse_args()
        main(args)
