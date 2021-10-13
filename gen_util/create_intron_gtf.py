# Takes a intron gtf file and corresponding appris txt file
# http://appris-tools.org/#/downloads
# Only returns appris principal 1 transcript(s)
# Generate intron gtf file by following these instructions
# https://groups.google.com/a/soe.ucsc.edu/g/genome/c/VYBl3k3IX4I

from collections import defaultdict
import argparse
import gzip


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Generates a intron gtf or translates ENSMBL transcript IDs to genesym')
	
	parser.add_argument("command",
			metavar="<command>",
			help="Please choose 'gtf' or 'genesym'")
			
	parser.add_argument('--appris', required=True,
			help='Appris .txt file. See http://appris-tools.org/#/downloads')
	
	parser.add_argument('--gtf', required=False,
			help='Intronic gtf file (used with <command> = "gtf"). See https://groups.google.com/a/soe.ucsc.edu/g/genome/c/VYBl3k3IX4I')
			
	parser.add_argument('--tsv', required=False,
			help='Tsv file (used with <command> = "genesym"). Make sure "ENST" ids in 1st column. Write to STDOUT')
	args = parser.parse_args()
	
	appris_file = args.appris
	
	if args.command == 'gtf':
		gtf_file = args.gtf
		# Import appris file
		# appris_dict[geneID][appris_tag] = [transcriptIDs with same appris_tag]
		appris_dict = defaultdict(lambda: defaultdict(lambda: []))

		with open(appris_file, 'r') as f:
			for line in f:
				line = line.strip().split('\t')
				transcriptID = line[2]
				geneID = line[1]
				tag = line[-1]
				appris_dict[geneID][tag].append(transcriptID)


		# Now decide on the principal isoform to use
		# Always report PRINCIPAL:1 if it exists
		# Else report PRINCIPAL:2, 3, 4, or 5
		# appris_p1[transcriptIDs] = geneID		
		appris_p1 = defaultdict(lambda: '')
		for geneID in appris_dict.keys():
			tag_list = appris_dict[geneID].keys()
			counter = 1
			tag = "PRINCIPAL:" + str(counter)
			# Figure out what is best tag
			while appris_dict[geneID][tag] == []:
				counter += 1
				tag = "PRINCIPAL:" + str(counter)
				if counter > 6:
					break
			for transcriptID in appris_dict[geneID][tag]:
				appris_p1[transcriptID] = geneID
				
		# Import GTF File 
		# Export only import appris_p1 transcripts
		# Some genes may have > 1 P1 transcripts. These transcripts may also share a significant amount of introns
		# Reduce the file down to only unique intronic positions
		# intron[(chrm,start,stop)] = boolian 0 (False) 1 (True) if the positional tuple has been seen before
		intron_pos = defaultdict(lambda: 0)

		if gtf_file.endswith('.gz'):
			f = gzip.open(gtf_file, 'r')
			outname = gtf_file.rstrip('gtf.gz') + 'apprisP1.gtf.gz'
		else:
			f = open(gtf_file, 'r')
			outname = gtf_file.rstrip('gtf') + 'apprisP1.gtf.gz'

		outfile = gzip.open(outname, 'wt')

		for line in f:
			line = line.strip().split('\t')
			attribute = line[-1].split('; ')
			# Noticed that intron gtf file has the transcriptID
			# for the geneID. Need to correct this
			transcriptID = attribute[1].split('"')[1].split('.')[0]
			geneID = appris_p1[transcriptID]
			
			# Store position of introns
			chrm = line[0]
			start = line[3]
			stop = line[4]
			position = (chrm, start, stop)
			
			if (geneID != '') and (intron_pos[position] == 0):
				#gene_attribute = 'gene_id "' + geneID + '"'
				#new_attribute = [gene_attribute] + attribute[1:]
				#new_attribute = '; '.join(new_attribute)
				#outfile.write('%s\t%s\n' % ('\t'.join(line[:-1]), new_attribute))
				outfile.write('%s\n' % ('\t'.join(line)))
				intron_pos[position] = 1
	
	elif args.command == 'genesym':
		# Import appris file
		# appris_dict[transcriptID] = geneID
		appris_dict = defaultdict(lambda: '')

		with open(appris_file, 'r') as f:
			for line in f:
				line = line.strip().split('\t')
				transcriptID = line[2]
				geneSym = line[0]
				appris_dict[transcriptID] = geneSym
				
		
		tsv_file = args.tsv
		
		with open(tsv_file, 'r') as f:
			for line in f:
				line = line.strip().split('\t')
				transcriptID = line[0]
				geneSym = appris_dict[transcriptID.split('.')[0]]
				
				if geneSym != '':
					new_text = line[0].replace(transcriptID.split('_')[0], geneSym)
					line[0] = new_text
					print('\t'.join(line))
				else:
					print('\t'.join(line))
				
				
				
					
	
