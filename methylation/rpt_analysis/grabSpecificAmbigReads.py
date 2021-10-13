### Grab reads with specific # of multimappin locations from *ambig.bam file
### returns sam file of reads

import argparse
import pysam
from collections import defaultdict

def main(args):

	
	# Parse arguments
	bam = args.bam
	limit = args.limit

	orig_reads = defaultdict(lambda: ())
	bad_set = set()

	# Parse thru bam
	# I ASSUME that every unique read has at least 2 multimapping locations 

	bamfile = pysam.AlignmentFile(bam,"rb")
	for line in bamfile:
		current_name = line.qname
		if current_name in bad_set:
			continue
		orig_reads[current_name] = orig_reads[current_name] + (line,)
		# This just helps with processing speed. If a read has too many locations, then it is dropped from the 
		# dictionary and no longer considered
		# In the above if statement, if that read name pops up again, ignore it altogether
		if len(orig_reads[current_name]) > limit:
			orig_reads.pop(current_name, None)
			bad_set.add(current_name)

	output = bam + ".SPECIFIC%sloc.bam" % limit
	outfile = pysam.AlignmentFile(output, 'wb', template=bamfile)

	for key in orig_reads.keys():
		if len(orig_reads[key]) == limit:
			for line in orig_reads[key]:
				outfile.write(line)
	outfile.close()
	bamfile.close()


if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='See comments in code')
        parser.add_argument('-b', dest='bam', type=str, help='*ambig_bam file from bismarkMultimappingLocations')
	parser.add_argument('-l', dest='limit', default = 1000, type=int, help='Optional: Limit on the # of multimapped locations a read can have (takes integer and run time is longer). If a read has more locations than the specified limit, it is NOT analyzed')
	args=parser.parse_args()
        main(args)
