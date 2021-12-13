# Takes a 10x bam file and a list of cellbarcodes to extract (can be gzipped)
# Ex: python split_10x_cellbarcode.py 10x.bam cellbarcodes.txt

import pysam
import gzip
import os
from sys import argv
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--bam","-b", type = str, dest = "bam", help = "Input cellranger count outputted bam file")
parser.add_argument("-c", "--barcodes", type = str, dest = "barcodes", help = "Only output reads associated with specified barcodes (txt or gzip format)")
parser.add_argument("--exon","-e", type = str, dest="exon", default="n", help = "Type 'y' if only want to output read pairs that overlap 50% with exons")
args = parser.parse_args()

# Read userbarcodes into dictionary
userbars = defaultdict(lambda: 0)

if args.barcodes.endswith('.gz'):
	with gzip.open(args.barcodes, 'r') as f:
		for line in f:
			line = line.strip()
			userbars[line] += 1
else:
	with open(args.barcodes, 'r') as f:
		for line in f:
			line = line.strip()
			userbars[line] += 1

# Parse thru bam file and store lines into dictionary
	# This can be quite memory intensive so be careful
	# Multithreading will also improve the speed of this

bamlines = defaultdict(lambda: [])

with pysam.AlignmentFile(args.bam, "rb") as f:
	tenxheader = f.header
	for n, line in enumerate(f):
	#	if n % 1000000 == 0:
	#		print(n)

		# This checks if line has valid cell barcode or not
		cbcheck = 0

		try:
			chrm = "chr" + str(line.reference_name)
			tags = line.get_tags()
			for t in tags:
				if t[0] == 'CB':
					CBcall = t[1]
					cbcheck += 1
					break
			if cbcheck == 0:
				continue

			if args.exon == 'y':
				for t in tags:
					if t[0] == 'RE' and t[1] == 'E':
						cbcheck += 1
						break
				if cbcheck != 2:
					continue			

			outstring = line.to_string()
			mapq = line.mapping_quality
			flag = int(line.flag)

			# Filter by cellbarcodes first
			if userbars[CBcall] >= 1:
				# Filter by mapping quality
				#if int(mapq) == 255:
				# Filter by flag
					#if flag not in [77, 141, 339, 355, 403, 419]:
				bamlines[CBcall].append(line)
		except:
			pass

print('Writing to output')
# Output into sams without headers :(
for bc in bamlines.keys():
	os.mkdir('./' + bc)
	outname = './' + bc + '/align.bam'
	with pysam.AlignmentFile(outname, "wb", header=tenxheader) as outf:	
		for i in bamlines[bc]:
			outf.write(i)

	print(bc)
