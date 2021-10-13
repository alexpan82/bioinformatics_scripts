# Turns bam file into a bed file

import pysam
from sys import argv

script, bam = argv

outname = bam.rstrip('bam')
outname = outname + 'bed'
outfile = open(outname, 'w')
outfile.write('#chr\tstart\tstop\tname\n')
with pysam.AlignmentFile(bam, "rb") as f:
	for line in f:
		start = str(line.pos)
		stop = str(int(start) + int(line.query_alignment_length) + 1)
		chrm = line.reference_name
		name = line.qname
		outfile.write('%s\t%s\t%s\t%s\n' % (chrm, start, stop, name))
	f.close()

outfile.close()
