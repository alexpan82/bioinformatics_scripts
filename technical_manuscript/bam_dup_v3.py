import pysam
import argparse
import pybedtools
from collections import defaultdict


def main(args):


	bam = args.bam
	option = args.option

	if option == 'pos':

		#pos = set()
		# pos[bam_name] = set("R1 chrm start stop", "R2 chrm start stop")
		pos = defaultdict(lambda: [])

		# Read bam file and store mates in pos dictionary
		#print 'Reading thru bam'
		with pysam.AlignmentFile(bam, "rb") as f:
			for n, line in enumerate(f):
				name = line.qname
				start = int(line.pos)
				stop = int(start) + len(line.query_sequence)
				chrm = line.reference_name
				pos[name].append([chrm, start, stop])

		#uniqset = set()

		# Add "chr\tstart\tstop,chr\tstart\stop" to set
		# Make sure to sort the string by start position
		# Also calculate avg frag length by storing in list
		
		fraglen = []
		#countdup = defaultdict(lambda: 0)
		readnamepos = defaultdict(lambda: [])
		for name in pos.keys():
			if pos[name][0][1] >= pos[name][1][1]:
				index = 0
			elif pos[name][0][1] < pos[name][1][1]:
				index = 1
			if index == 0:
				addstring = '%s,%s' % ('\t'.join(map(str,pos[name][1])), '\t'.join(map(str, pos[name][0])))
				fraglen.append(pos[name][0][2]-pos[name][1][1])
			elif index == 1:
				addstring = '%s,%s' % ('\t'.join(map(str,pos[name][0])), '\t'.join(map(str, pos[name][1])))
				fraglen.append(pos[name][1][2]-pos[name][0][1])

			#uniqset.add(addstring)
			#countdup[addstring] += 1
			readnamepos[addstring].append(name)
		

		fragavg = sum(fraglen) / len(fraglen)
		print("%s avg len: %s" % (bam, fragavg))

		n += 1
		n /= 2
		#num_uniq = len(uniqset)
		num_uniq = len(readnamepos.keys())
		#print('%s\t%s\t%s\t%s' % (bam, num_uniq, n, (float(n)-num_uniq)/float(n)))

		#outdup = open('%s.dup.pos.bed' % '_'.join(bam.split('_')[0:4]), 'w')
		#outuniq = open('%s.uniq.pos.bed' % '_'.join(bam.split('_')[0:4]), 'w')
		
		numreadsdup = 0
		numreadsuniq = 0

		for i in readnamepos.keys():
			if len(readnamepos[i]) > 1:
				for name in readnamepos[i]:
					#outdup.write('%s\t%s\n%s\t%s\n' % (i.split(',')[0], name, i.split(',')[1], name))
					numreadsdup += 1
				#readnamepos[i] = readnamepos[i][1:]
			else:
				for name in readnamepos[i]:
					#outuniq.write('%s\t%s\n%s\t%s\n' % (i.split(',')[0], name, i.split(',')[1], name))
					numreadsuniq += 1
				#readnamepos[i] = readnamepos[i][1:]
		#print('Number of reads NOT position duplicated:\t%s' % numreadsuniq)
		#print('Number of reads position duplicated:\t%s' % numreadsdup)

	elif option == 'seq':

		#seq = set()
		# seq[bam_name] = set([[R1chrm, start, stop, line.query_sequence], [R2chrm, start, stop, line.query_sequence]])
		seq = defaultdict(lambda: [])

		# Assumes paired reads in bam file have the same name and R1 comes before R2 in the bam file
		#print 'Reading thru bam'
		with pysam.AlignmentFile(bam, "rb") as f:
			bed = ''
			for n, line in enumerate(f):
				name = line.qname
				start = int(line.pos)
				stop = int(start) + len(line.query_sequence)
				chrm = line.reference_name
				seq[name].append([chrm, start, stop, line.query_sequence])

		uniqset = set()

		for name in seq.keys():
			if seq[name][0][1] >= seq[name][1][1]:
				index = 0
			elif seq[name][0][1] < seq[name][1][1]:
				index = 1
			if index == 0:
				addstring = '%s,%s' % (seq[name][1][3], seq[name][0][3])
			elif index == 1:
				addstring = '%s,%s' % (seq[name][0][3], seq[name][1][3])

			uniqset.add(addstring)

	
		n += 1
		n /= 2
		num_uniq = len(uniqset)
		print('%s\t%s\t%s\t%s' % (bam, num_uniq, n, (float(n)-num_uniq)/float(n)))



	else:
		print('Please type -o pos or -o seq')
		quit()


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='See comments in code')
	parser.add_argument('-b', dest='bam', type=str, help='Uniquely aligned unsorted PE bam')
	parser.add_argument('-o', dest='option', type=str, help='Type "pos" or "seq" to calculate positional or sequence duplication, respectively')
	args=parser.parse_args()
	main(args)






