# This script will pull out the 3' sequence AFTER the end of the J gene (if it exists)
# This script has 2 modes: pull and paste
# Run like this:
	# 1) python igblastClass.py pull [CHANGEO TSV] > [FASTA]
	# 2) blastn -db nt -outfmt "6 qseqid bitscore stitle" -query [FASTA] -out [OUT]
	# 3) python igblastClass.py paste [CHANGEO TSV] [OUT] > [FINAL]

def imgtToVDJ(string):
	string = string.split('.')
	return ''.join(string)

from sys import argv
from collections import defaultdict

# Read thru tsv and output 3' end after J-gene
def pull3prime(changeoTsv):

	with open(changeoTsv, 'r') as f:
		for line in f:
			try:
				line = line.strip().split('\t')
				NAME = line[0]
				SEQUENCE = line[1]

				# Where J ends relative to query sequence
				J_SEQ_END = int(line[24]) + int(line[25]) - 1

				# print only 3' end after J
				# Dont print if too short

				threeprime_keep = SEQUENCE[J_SEQ_END:]
				if len(threeprime_keep) > 5:
					print('>%s\n%s' % (NAME, threeprime_keep))


			except:
				pass

# Strings to search for in blast output. Modify however necessary
def possibleIG(character):
	a = 'ig' + character
	b = 'igh' + character
	c = 'immunoglobulin heavy constant ' + character
	d = 'immunoglobulin ' + character
	e = 'ig ' + character
	return [a, b, c, d, e]

# Associate CHANGEO IGBLAST lines with class
def pasteTsv(changeoTsv, BLAST):

	querydic = defaultdict(lambda: "")
	maxbitscore = defaultdict(lambda: 0)

	# Read thru blast output and find strings. Only read from highest bit value
	with open(BLAST, 'r') as f:
		for line in f:
			line = line.strip().split('\t')
			try:
				bitscore = int(line[1])
				query = line[0]
				alignment = line[2].lower()

				if bitscore >= maxbitscore[query]:

					maxbitscore[query] = bitscore

					for c in ['a', 'd', 'm', 'g', 'e']:
						for search in possibleIG(c):
							if search in alignment:
								querydic[query] = possibleIG(c)[0].upper()
								break
							else:
								pass
			except:
				pass

	# Open changeo file and append class to last column
	with open(changeoTsv, 'r') as f:
		for n, line in enumerate(f):
			if n == 0:
				print(line.strip() + "\tIG_CLASS")
			else:
				current_query = line.strip().split('\t')[0]
				
				if querydic[current_query] == '':
					current_query = current_query.split()[0]


				print(line.strip() + '\t%s' % (querydic[current_query]))

# Main

if argv[1] == 'pull':
	try:
		pull3prime(argv[2])
	except:
		print("Please input a CHANGEO (MakeDb) formatted table from igBLAST output as the 2nd argument")
		quit()
elif argv[1] == 'paste':
	try:
		pasteTsv(argv[2], argv[3])
	except:
		print("""Please input:
			2nd arg: CHANGEO (MakeDb) formatted table from igBLAST output as the 2nd argument
			3rd arg: 'nblast -db nt -outfmt "6 qseqid bitscore stitle"' output generated using the fasta file from igblastClass.py pull""")
		quit()
else:
	print("Please input 'pull' or 'paste' as the first argument")
	quit()
