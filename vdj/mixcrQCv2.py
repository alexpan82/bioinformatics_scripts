# Counts the number of aligned reads from Mixcr
# Takes Mixcr job logs as input

from sys import argv
from collections import defaultdict
#from scipy import stats

# Array to hold alignment info
alignments = defaultdict(lambda: defaultdict(lambda: 'NA'))

headers = ['Total sequencing reads:', 'Successfully aligned reads:', 
'Overlapped:', 'Overlapped and aligned:', 
'IGH chains:', 'IGK chains:', 'IGL chains:', 
'TRA chains:', 'TRB chains:', 'TRD chains:', 'TRG chains:', 'TRD chains:', 'TRA,TRD chains:',
'Alignments already with CDR3 (no overlapping is performed):',
'Final clonotype count:',
'Reads used in clonotypes, percent of total:', 'Reads used in clonotypes before clustering, percent of total:' ,
'Reads dropped due to the lack of a clone sequence:', 'Reads dropped due to low quality:',
'Reads dropped due to failed mapping:', 'Reads dropped with low quality clones:']

for infile in argv[1:]:
	with open(infile, 'r') as f:
		for line in f:
			if any(h in line for h in headers):
				line = line.strip().split(':')
				#print(line)
				if alignments[infile]['Alignments already with CDR3 (no overlapping is performed)'] != 'NA' and line[0].strip() == 'Alignments already with CDR3 (no overlapping is performed)':
					alignments[infile]['Alignments already with CDR3 (no overlapping is performed)'] += ', %s' % line[1].strip()
				else:
					alignments[infile][line[0].strip()] = line[1].strip()

			if 'Reads dropped with low quality clones' in line[0]:
				break



# Create tsv
outfile = open('mixcrQCTable.tsv', 'w')
hstring = '\t'.join(list(map(lambda x: x.rstrip(':'), headers)))
outfile.write('\t%s' % hstring)

for sample in sorted(alignments.keys()):
	outfile.write('\n%s' % (sample))
	for h in headers:
		h = h.rstrip(':')
		outfile.write('\t%s' % alignments[sample][h])
	
outfile.close()
