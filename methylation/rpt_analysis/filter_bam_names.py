### given a list of bam read names, outputs all reads that are NOT in that list
### returns sam file of reads

import argparse
import pysam
from collections import defaultdict

def main(args):

	
	# Parse arguments
	bam = args.bam
	blist = args.blist

	bad_set = set()

	# Read list
	with open(blist, 'r') as f:
		for line in f:
			bad_set.add(line.strip())

	output = bam.strip('bam')
	output = output + "RPTFILTEREDOUT.sam"
	outfile = open(output, 'w')
	# Parse thru bam
	with pysam.AlignmentFile(bam, "rb") as f:
		for line_num,line in enumerate(f):
			current_name = line.qname
			if current_name not in bad_set:
				outfile.write("%s\n" % line.to_string())
			else:
				continue

	f.close()
	outfile.close()


if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='See comments in code')
        parser.add_argument('-b', dest='bam', type=str, help='*ambig_bam file from bismarkMultimappingLocations')
	parser.add_argument('-l', dest='blist', type=str, help='list of bam read names')
	args=parser.parse_args()
        main(args)
