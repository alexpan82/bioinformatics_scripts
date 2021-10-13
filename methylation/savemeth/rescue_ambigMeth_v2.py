### Written by: Alex Pan
### Rescue methylation calls from multimapped PE reads
### Multimapped PE reads are aligned in SE mode
### Returns uniquely mapped SE reads
### There are 2 types of reads this script can save
	### 1) Both mates uniquely map close to each other (SE alignment)
	### 2) R1 (or R2) uniquely maps, and R2 multimaps near the R1 SE alignment

### File requirements
	### 1) Uniquely aligned SE bam (generated from multimappe PE reads)
	### 2) Fastq file used to generate (1)
	### 3) *ambig.bam file containing PE multimapping locations (see BismarkMultimappingLocations)
	### 4) Reference genome .fa file

# Dependancies:
	# pysam
	# pybedtools

import pysam
import argparse
import pybedtools
import gzip
from collections import defaultdict
import string
import sys
import re

################################################################################################################################################***
# Reverse complement 
def revcomp(seq):
	trans = string.maketrans('ATGC','TACG')
	return seq.translate(trans)[::-1]

################################################################################################################################################***
# Complement
def comp(seq):
	trans = string.maketrans('ATGC','TACG')
	return seq.translate(trans)

################################################################################################################################################***
# Progress bar lol
def progress(count, total, status=''):
	bar_len = 60
	filled_len = int(round(bar_len * count / float(total)))

	percents = round(100.0 * count / float(total), 1)
	bar = '=' * filled_len + '-' * (bar_len - filled_len)

	sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
	sys.stdout.flush()

################################################################################################################################################***
# Get line count of file for purpose of progress bar
def file_len(fname):
	print 'Counting lines in %s...' % fname
	if fname.endswith('bam'):
		with pysam.AlignmentFile(fname, "rb") as f:
			for i, l in enumerate(f):
				pass
			f.close()
	elif fname.endswith('gz'):
		with gzip.open(fname, "r") as f:
			for i, l in enumerate(f):
				pass
			f.close()
	else:
		with open(fname, "r") as f:
			for i, l in enumerate(f):
				pass
			f.close()
	return i + 1

################################################################################################################################################***
# Function to deal with insertions and deletions in reads. A key part of this code is translating the relative position of CpGs on a reference sequence to the position on the
# bisulfite converted read. However the read can have indels which change the expected position(s). This function returns a value that reflects
# how much these positions have changed due to indels. Needs cigar tuple and the reference CpG relative postion (index)
# Cigar flags are formatted as follows: ((flag, position), ...)
	# 0 denotes normal sequence from last position to current position
	# 1 denotes insertion(s) as described by the position
	# 2 denotes deletion
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
	#index -= 1
	if index < 0:
		print seq
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

# Call methylation from methylation flag. 'Z' is methylated, 'z' is unmethylated
# Takes a sequence, index and itr (see above)
def callmethflag(seq, index):
	#print index
	#print seq
	if seq[index] == 'Z':
		meth = 1
	elif seq[index] == 'z':
		meth = -1
	else:
		meth = 'wut'

	return meth

################################################################################################################################################***

def main(args):


	bam = args.bam
	ambigbam = args.ambigbam
	fastq = args.fastq
	fasta = args.fasta
	outname = args.outname
	
	################################################# Parse thru files, set up appropriate data structures, group reads by how they can be saved #################################################

	# First parse thru unique SE bam file to obtain reads, sequence information, cigar string, and methylation call string
	uniqbamdic = defaultdict(lambda: ())

	# Holds names of R1 and R2 reads
	uniqr1 = set()
	uniqr2 = set()
	
	print "Parsing through %s ..." % bam

	with pysam.AlignmentFile(bam, "rb") as f:
		
		for line in f:
			start = int(line.pos)
			seq = line.query_sequence
			stop = int(start) + len(seq)
			chrm = line.reference_name
			current_name = line.qname
			cigar = line.cigartuples
			call = line.get_tags()

			# only want alignments to chr
			if chrm.startswith('chr'):
				uniqbamdic[current_name] = (chrm, start, stop, seq, cigar, call)

				# Add names to these sets to figure out which reads we need to look for in ambigbam
				if '_1:N:' in current_name:
					uniqr1.add(current_name)
				elif '_2:N:' in current_name:
					uniqr2.add(current_name)
		f.close()

	#############################################################################

	# Now figure out which reads have mate uniquely mapped in close proximity (we're working with PE 150 data, so if the mate is 1000bp away, I count it)


	r2tmp = set(map(lambda x: x.replace('_2:N:', '_1:N:', 1), list(uniqr2)))
	
	# set of R1 reads whos' mates also uniquely mapped
	r1mates = uniqr1 & r2tmp

	del r2tmp

	# set of R2 reads whos' mates also uniquely mapped (Same as r1mates)
	r2mates = set(map(lambda x: x.replace('_1:N:', '_2:N:', 1), list(r1mates)))

	# Loop thru r1mates checking which R2 mates map within 1000bp
	# Store these reads in r1proxmates and r2proxmates
	r1proxmates = set()
	r2proxmates = set()
	
	count = 0
	for r1name in r1mates:
		r2name = r1name.replace('_1:N:', '_2:N:', 1)
		chr1 = uniqbamdic[r1name][0]
		chr2 = uniqbamdic[r2name][0]
		start1 = uniqbamdic[r1name][1]
		start2 =uniqbamdic[r2name][1]
		if chr1 == chr2:
			if (start2-1000) <= start1 and start1 < (start2+1000):
				r1proxmates.add(r1name)
				r2proxmates.add(r2name)

	# Uniquely mapped R1 reads without R2 mate and uniquely mapped R2 reads without R1 mate 
	# Will have to search thru *ambig.bam to find where those reads multimapped
	r1womate = uniqr1-r1mates
	r2womate = uniqr2-r2mates

	#############################################################################

	# Search thru ambigbam for mates that multimapped
	# Store read location info 

	print "Parsing through %s for reads that multimapped near uniquely aligned mate ..." % ambigbam

	# Count some things
	matefar = set()
	allambig = set()

	multimapmates = defaultdict(lambda: ())
	#total = file_len(ambigbam)
	# Get line count of bam file
	
	with pysam.AlignmentFile(ambigbam, "rb") as f:
		bed = ''
		for n,line in enumerate(f):
			#progress(n, total, status='Progress')
			mstart = int(line.pos)
			mseq = line.query_sequence
			mstop = int(mstart) + len(mseq)
			mchr = line.reference_name
			multimapped_name = line.qname
			allambig.add(multimapped_name)
			cigar = line.cigartuples
			if mchr.startswith('chr'):
				if '_1:N:' in multimapped_name:
					uniqmapped_mate = multimapped_name.replace('_1:N:', '_2:N:', 1)
					if uniqmapped_mate in r2womate:
						matefar.add(multimapped_name)
						uchr = uniqbamdic[uniqmapped_mate][0]
						ustart = uniqbamdic[uniqmapped_mate][1]
						if mchr == uchr:
							if (mstart-1000) <= ustart and ustart < (mstart+1000):
								multimapmates[multimapped_name] = (mchr, mstart, mstop, cigar)
								# Continuously add to string
								bed = bed + '%s\t%s\t%s\t%s\n' % (mchr, mstart, mstop, multimapped_name)

				# I have these in if statements so script doesnt have to check r1womate AND r2womate for string
				elif '_2:N:' in multimapped_name:
					uniqmapped_mate = multimapped_name.replace('_2:N:', '_1:N:', 1)
					if uniqmapped_mate in r1womate:
						matefar.add(multimapped_name)
						uchr = uniqbamdic[uniqmapped_mate][0]
						ustart = uniqbamdic[uniqmapped_mate][1]
						if mchr == uchr:
							if (mstart-1000) <= ustart and ustart < (mstart+1000):
								multimapmates[multimapped_name] = (mchr, mstart, mstop, cigar)
								# Continuously add to string
								bed = bed + '%s\t%s\t%s\t%s\n' % (mchr, mstart, mstop, multimapped_name)
			else:
				pass
		f.close()

	print "Number of uniquely mapped reads who's mate also uniquely mapped:\t%s" % (len(r1mates)*2)
	print 'Mate multimapped far away from uniquely mapped read:\t%s' % len(matefar-set(multimapmates.keys()))
	print 'Mate multimapped close from uniquely mapped read:\t%s' % len(multimapmates.keys())
	print 'Mate not uniquely mapped or multimapped:\t%s' % len(set(r1womate | r2womate) - allambig)
	del (matefar, allambig)

	#############################################################################
	# Need to grab the reference sequence for the location of the multimapped read so that we may call methylation
	# THis is slow step. make separate loop to store values in a dic

	genseqdic = {}
	print "Running bedtools getfasta on passed multimapped reads' locations ...\n"
	pybedobject = pybedtools.BedTool(bed, from_string = True)
	pybedobject = pybedobject.sequence(fi=fasta, name=True)
	del bed
	for line in open(pybedobject.seqfn):
		line = line.strip()
		if line.startswith('>') is True:
			readname = line.lstrip('>')
		else:
			genseqdic[readname] = line.upper()
	del pybedobject

	# Now that we have which multimapped mates multimap in the vincinity of the uniquely mapped mate, we need to grab sequence info from the fastq file
	# This is bcuz *ambig.bam erases methylation info from the sequence
	# nultimapmates[name] = (chr, start, stop, fastq seq)
	tmpset = set(multimapmates.keys())

	print "Obtaining fastq sequence info for passed multimapped reads ..." 
	#total = file_len(fastq)
	i = 1
	with gzip.open(fastq, 'rb') as f:
		for n,line in enumerate(f):
			#progress(n, total, status='Progress')
			line = line.strip()
			if i == 1:
				name = line.lstrip('@')
			elif i == 2:
				multimapmates[name] = multimapmates[name] + (line,)
			elif i ==4:
				i = 0
			i += 1

		f.close()

	for i in set(set(multimapmates.keys())-tmpset):
		multimapmates.pop(i)

	del tmpset
	print '\nRecovering methylation ...'

	# Parse optional arguments... Here, the names of reads that have initially PE ambiguously or un- aligned are given so that they may be differentiated
		# From here, we will store names and cpgs postions once the time comes

	PEanames = args.anames
	PEunnames = args.unnames
	PEnonuniq_dic = {}
	if PEanames != 'na':
		print 'Additionally writing output for CpGs from specified reads from %s only' % PEanames
		# set to hold names
		with open(PEanames, 'r') as f:
			for line in f:
				line = line.strip()
				PEnonuniq_dic[line] = 'pe_ambig'
			f.close()

	if PEunnames != 'na':
		print 'Additionally writing output for CpGs from specified reads from %s only' % PEunnames
		PEunnames_dic = {}
		with open(PEunnames, 'r') as f:
			for line in f:
				line = line.strip()
				PEnonuniq_dic[line] = 'pe_unalign'
			f.close()


	############################################################################## Call methylation on reads ##############################################################################

	# Let's start with reads where only 1 mate is uniquely mapped, but other mate has multimapped location in vicinity
	# The uniquely mapped read is in the orientation of the genome 5' to 3'
	
	# Looking at the alignment flag for the uniquely mapped mate will tell us about how the multimapped mate should align
	# The mate may also be read as the reverse comp to the orientation of the genome, so need to be able to call CpGs accordingly

	# Open a temp file to contain methylation information
	tmpoutput = outname + ".tmp"
	tmpopen = open(tmpoutput, 'w')

	for mname in multimapmates.keys():
		try:
			if '_1:N:' in mname:
				uname = mname.replace('_1:N:', '_2:N:', 1)
			elif '_2:N:' in mname:
				uname = mname.replace('_2:N:', '_1:N:', 1)

			# Get mapping info from uniquely mapped mate
				# xr is the conversion flag for the sequence, xg is the converted genome xr aligned to
				# xr=ct and xg=ct or xr=ga and xg=ga means read was 5' to 3' in fastq
					# Gives info on top strand
				# xr=ga and xg=ct or xr=ct and xg=ga means read was 3' to 5' in fastq
					# Gives info on bottom strand
			# In fastq files, the orientation of R2 is the reverse comp of R1
			u_xr = uniqbamdic[uname][5][3][1]
			u_xg = uniqbamdic[uname][5][4][1]

			# Bam sequences and methylation call flags are always reported top strand (5' to 3') of the reference genome, no matter which
				# strand the bisulfite coverted read may have come from
			# Orient multimapped sequence to the top strand (5' to 3') of the reference genome so that we may normalize seqs of our multimapped reads
				# to reported bam seq and methylation calls of our uniquely aligned reads
			# Fastq seq are always 5' to 3' (R2 is just on opposite strand so it's orientation is the reverse complement or R1)
			# Depending on which bisulfite converted strand the fastq seq is, it may not be in this top strand (5' to 3') orientation
			# Alignments to OT and CTOT yield meth info on top strand. OB and CTOB yield info about bottom strand

			mseq = multimapmates[mname][4]
		
			if u_xr == 'CT' and u_xg == 'CT': # OT
				# The uniquely aligned read, in it's original fastq orientation, was 5' to 3' and mapped to the 5' to 3' C to T genome
				# Since the fastq of the mate is in the revcomp orientation, just need to take revcomp
				mseq = revcomp(mseq)
				strandinfo = '+'
			elif u_xr == 'GA' and u_xg == 'GA': # CTOB
				mseq = revcomp(mseq)
				strandinfo = '-'
			elif u_xr == 'CT' and u_xg == 'GA': # OB
				# In these cases, CTOT or OB is what aligned. The bam file will contain the complement of the fastq seq, as the fastq seq of CTOT or OBSince bams only report seq 5' to 3' of top strand,
				# the aligned mate is reported as the revcomp what was in the fastq. Therefore, the multimapped mate must have been
				# originally reported in the top strand orientation in the fastq, so no modification is necessary
				strandinfo = '-'

			elif u_xr == 'GA' and u_xg == 'CT': # CTOT
				strandinfo = '+'

			# Find locations on multimapped fastq (that is now in 5' to 3' top) that are in CG context
			# Do this by finding CGs on reference seq
			rel_cg_loc_ref = sorted(set([m.start() for m in re.finditer('CG', genseqdic[mname])]))

			# The issue now is that indels in the fastq sequence will shift these relative positions

			for loc in rel_cg_loc_ref:
				# Relative positions start with zero but genome pos start at 1
				shift = indel(multimapmates[mname][3], loc)
				meth = callmeth(mseq, loc, shift)
				mchr = multimapmates[mname][0]
				gloc = loc + multimapmates[mname][1] + 1
				if strandinfo == '-':
					gloc += 1
				# Write chr, location, methylation (1 || -1), strand (+ || -), read name, origin of read (pe_ambig || pe_unalign), recovery method (se_ambig, se_align)
				tmpopen.write('%s\t%s\t%s\t%s\t%s\t%s\tse_ambig\n' % (mchr, gloc, meth, strandinfo, mname, PEnonuniq_dic[mname]))
		except:
			print 'Wut\t%s' % mname

	# Now let's call methylation on uniquely aligned reads. There are 2 subsets
		# 1) Reads who uniquely aligned but mate did not map at all
		# 2) Read and mate both uniquely aligned
	# Names of the uniquely aligned reads are in uniqr1 and uniqr2
	uniqreads = uniqr1 | uniqr2
	# We dont need these guys anymore	
	del (uniqr1, uniqr2)

	# Recovery is easy as we just look at the methylation call string!
	
	for uname in uniqreads:
		methflags = uniqbamdic[uname][5][2][1]
		# Find out whether seq gives info about top or bot strandgdic = defaultdict(lambda: defaultdict(lambda: (0,0)))

		u_xr = uniqbamdic[uname][5][3][1]
		u_xg = uniqbamdic[uname][5][4][1]
		if u_xr == 'CT' and u_xg == 'CT': # OT
			strandinfo = '+'
		elif u_xr == 'GA' and u_xg == 'GA': # CTOB
			strandinfo = '-'
		elif u_xr == 'CT' and u_xg == 'GA': # OB
			strandinfo = '-'
		elif u_xr == 'GA' and u_xg == 'CT': # CTOT
			strandinfo = '+'

		# Find position of all CpG sites (z or Z)
		rel_cg_loc_ref = sorted(set([m.start() for m in re.finditer('Z', methflags, re.IGNORECASE)]))
		
		# Meth flags already account for indels
		for loc in rel_cg_loc_ref:

			meth = callmethflag(methflags, loc)
			mchr = uniqbamdic[uname][0]
			gloc = loc + uniqbamdic[uname][1] + 1
			if strandinfo == '-':
				gloc += 1
			tmpopen.write('%s\t%s\t%s\t%s\t%s\t%s\tse_align\n' % (mchr, gloc, meth, strandinfo, uname, PEnonuniq_dic[uname]))
	
	tmpopen.close()
	del PEnonuniq_dic
	print "Parsing through temporary file ..."

	# Combine files and write final output

	outfile = open(outname, 'w')

	# A CpG that both mates overlap should only be covered once
		# Or else we double-count coverage
	# So we need to keep track of what CpGs are covered by what read
	# In addition, both mates' methylation call for a CpG must agree. Else we throw away that CpG

	# cpgdic[readname]['chr\tcpg_position'] = (meth, strand)

	cpgdic = defaultdict(lambda: defaultdict(lambda: (0,0)))

	# Counts # of CpGs where mate pairs disagree on call
	disagree = 0
	disagree_list = defaultdict(lambda: 0)

	with open(tmpoutput, 'r') as f:
		for n,line in enumerate(f):
			line = line.strip().split('\t')
			try:
				# Skip lines with mutations
				if int(line[2]) == 1 or int(line[2]) == -1:
					generalname = line[4].split('_')[0]
					meth = int(line[2])
					chrm = line[0]
					strand = line[3]
					pos = line[1]
					loc = "%s\t%s" % (chrm, pos)
					# This loop checks whether the mates agree on the methylation call
					# If they don't, it's a sequencing error and we will ignore the methylation call entirely
					if cpgdic[generalname][loc] != (0,0):
						if cpgdic[generalname][loc] == (meth, strand):
							pass
						else:
							cpgdic[generalname].pop(loc)
							disagree += 1
							disagree_list[loc+'-'+generalname] = 1
					else:
						cpgdic[generalname][loc] = (meth, strand)
			except:
				pass
			#if n % 1000000 == 0:
				#print n

		f.close()

	print "Number of CpGs where mate pairs disagree on call:\t%s" % disagree

	# Since loc keys are a set, we have eliminated any overlaps between R1 and R2
	# Manipulate the dic such that it gives us the total methylation information

	print "Writing final output file ..."

	# Output dic struc: finalout[location] = [coverage, numCs]
	finalout = defaultdict(lambda: [0,0])
	count = set()
	for name in cpgdic.keys():
		for loc in cpgdic[name].keys():
			count.add(loc)
			finalout[loc][0] += 1
			if cpgdic[name][loc][0] == 1:
				finalout[loc][1] += 1

	print 'Unique CpG locations saved:\t%s' % len(count)

	for l in sorted(finalout.keys()):
		outfile.write("%s\t%s\t%s\n" % (l, '\t'.join(map(str, finalout[l])), (float(finalout[l][1]) / float(finalout[l][0]))))

	outfile.close()

	# Now print out some statistics
		# Find average coverage and pooled methylation for the final output
	totalcov = 0.0
	totalC = 0.0
	totalpos = 0.0
	for l in finalout.keys():
		totalpos += 1.0
		totalcov += float(finalout[l][0])
		totalC += float(finalout[l][1])
	print 'Stats for final output %s: ' % outname
	print 'Unique CpGs saved: %s' % len(finalout.keys())
	print '\tAverage Coverage: %s' % (totalcov / totalpos)
	print '\tAverage Pooled Methylation: %s' % (totalC / totalcov)
		
	
	# Print output for each type of CpG separately
	if PEanames != 'na' and PEunnames != 'na':
		del cpgdic
		print 'Additionally writing output for CpGs from specified reads from %s and %s' % (PEanames, PEunnames)
		# cpgdic[origin of read, recovery method][readname]['chr\tcpg_position'] = (meth, strand)
		cpgdic = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: (0,0))))
		with open(tmpoutput, 'r') as f:
			for n,line in enumerate(f):
				line = line.strip().split('\t')
				# Skip lines with mutations
				try:
					if int(line[2]) == 1 or int(line[2]) == -1:
						generalname = line[4].split('_')[0]
						meth = int(line[2])
						chrm = line[0]
						strand = line[3]
						pos = line[1]
						origin = line[5]
						rtype = line[6]	
						loc = "%s\t%s" % (chrm, pos)
						if disagree_list[loc+'-'+generalname] == 1:
							pass
						else:
							key = origin + ' ' + rtype
							cpgdic[key][generalname][loc] = (meth, strand)
				except:
					pass

				#if n % 100000 == 0:
					#print n
			f.close()
		#print 'hi1'				
		# Output dic struc: finalout[location] = [coverage, numCs]
		for category in cpgdic.keys():
			del finalout
			finalout = defaultdict(lambda: [0,0])
			count = set()
			for name in cpgdic[category].keys():
				for loc in cpgdic[category][name].keys():
					count.add(loc)
					finalout[loc][0] += 1
					if cpgdic[category][name][loc][0] == 1:
						finalout[loc][1] += 1

			print 'Unique CpG locations saved for %s:\t%s' % (category, len(count))
			#outname = outname.rstrip('permeth')
			#outname = outname + '%s.permeth' % category
			#outfile = open(outname+''.txt, 'w')
			#for l in sorted(finalout.keys()):
			#	outfile.write("%s\t%s\t%s\n" % (l, '\t'.join(map(str, finalout[l])), (float(finalout[l][1]) / float(finalout[l][0]))))

			#outfile.close()

		# Now print out some statistics
			# Find average coverage and pooled methylation for the final output
			totalcov = 0.0
			totalC = 0.0
			totalpos = 0.0
			for l in finalout.keys():
				totalpos += 1.0
				totalcov += float(finalout[l][0])
				totalC += float(finalout[l][1])
			print 'Stats for final output %s: ' % category
			print 'Unique CpGs saved: %s' % len(finalout.keys())
			print '\tAverage Coverage: %s' % (totalcov / totalpos)
			print '\tAverage Pooled Methylation: %s' % (totalC / totalcov)




	print "Done characterzing methylation on each unique CpG \n\nWorkflow complete!"


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='See comments in code')
	parser.add_argument('-b', dest='bam', type=str, help='Uniquely aligned SE bam (generated from multimappe PE reads)')
	parser.add_argument('-f', dest='fasta', type=str, help='Reference fasta name')
	parser.add_argument('-a', dest='ambigbam', type=str, help='*ambig_bam file from bismarkMultimappingLocations')
	parser.add_argument('-q', dest='fastq', type=str, help='Reference fastq name')
	parser.add_argument('-o', dest='outname', type=str, help='Output file name')
	parser.add_argument('-anames', dest='anames', type=str, default='na', help='Optional argument: Additionally output coverage and methylation values for CpGs only from PE ambiguously/multimapped reads. Takes text file with read names on each new line.')
	parser.add_argument('-unnames', dest='unnames', type=str, default='na', help='Optional argument: Additionally output coverage and methylation values for CpGs only from PE unaligned reads. Takes text file with read names on each new line.')
	args=parser.parse_args()
	main(args)
