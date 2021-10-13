### Grab every read EXCEPT those with specific # of multimappin locations from *ambig.bam file
### returns sam file of reads

import argparse
import pysam
from collections import defaultdict

def main(args):

	
	# Parse arguments
	bam = args.bam
	limit = args.limit

	orig_reads = defaultdict(lambda: ())

	# Parse thru bam
	# I ASSUME that every unique read has at least 2 multimapping locations 

	bamfile = pysam.AlignmentFile(bam,"rb")
	for line in bamfile:
		current_name = line.qname
		orig_reads[current_name] = orig_reads[current_name] + (line,)


	output = bam + ".EXCEPT%sloc.bam" % limit
	outfile = pysam.AlignmentFile(output, 'wb', template=bamfile)

	for key in orig_reads.keys():
		if len(orig_reads[key]) != limit:
			for line in orig_reads[key]:
				outfile.write(line)
	outfile.close()
	bamfile.close()


if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='See comments in code')
        parser.add_argument('-b', dest='bam', type=str, help='*ambig_bam file from bismarkMultimappingLocations')
	parser.add_argument('-l', dest='limit', default = 2, type=int, help='Grab every read EXCEPT those with specific # of multimappin locations from *ambig.bam file')
	args=parser.parse_args()
        main(args)
