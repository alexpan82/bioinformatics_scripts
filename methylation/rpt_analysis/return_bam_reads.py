### Given a list of read names, return reads and mates from bam file

import argparse
import pysam
from collections import defaultdict

def main(args):

	
	# Parse arguments
	bam = args.bam
	names = args.names

	orig_reads = defaultdict(lambda: ())
	names_set = set()


	with open(names, 'r') as f:
		for line in f:
			names_set.add(line.strip())
	f.close()
	# Parse thru bam 

	output = bam + ".SPECIFICNAMES.sam"
	outfile = open(output, 'w')

	with pysam.AlignmentFile(bam, "rb") as f:
		for line in f:
			current_name = line.qname
			if current_name in names_set:
				outfile.write(line.to_string())
	f.close()
	outfile.close()

if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='See comments in code')
        parser.add_argument('-b', dest='bam', type=str, help='*ambig_bam file from bismarkMultimappingLocations')
	parser.add_argument('-n', dest='names', default = 1000, type=str, help='List of bam read names')
	args=parser.parse_args()
        main(args)
