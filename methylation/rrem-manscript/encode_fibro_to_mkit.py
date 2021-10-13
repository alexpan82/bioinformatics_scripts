from sys import argv

script, encode = argv
print('chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT')
with open(encode, 'r') as f:
	for line in f:
		try:
			line = line.strip().split('\t')
			chr = line[0]
			pos = line[2]
			name = chr + '.' + pos
			if line[5] == '+':
				strand = 'F'
			else:
				strand = 'R'
			cov = line[-2]
			freqC = str(round(float(line[-1]), 2)) + "0"
			freqT = str(round(100.00-float(line[-1]), 2)) + "0"
			print('\t'.join([name, chr, pos, strand, cov, freqC, freqT]))
		except:
			pass
