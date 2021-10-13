import sys

def filldict(f):
	dictionary = {}
	header = f.readline()
	lines = f.readlines()
	for line in lines:
		cols = line.split()
		chrm = cols[1]
		pos = int(cols[2])
		cov = int(cols[4])
		freqC = float(cols[5]) / 100.0
		numC = int(round(cov * freqC))
		
		try:
			dictionary[chrm][pos] = str(cov) + '\t' + str(numC) + '\t'+ str(freqC)
		except:
			dictionary[chrm] = {}
			dictionary[chrm][pos] = str(cov) + '\t' + str(numC) + '\t'+ str(freqC)
		
	return dictionary

def writeoutput(out, filleddict):
	out.write("chr\tstart\tcov\tnumC\t5mC_ratio\n")
	for c in sorted(filleddict.keys()):
		for p in sorted(filleddict[c].keys()):
			out.write(str(c) + '\t' + str(p) +'\t' + str(filleddict[c][p]) + '\n')

	out.close()

def main():
	filename = sys.argv[1]
	outname = sys.argv[2]
	if filename is None or outname is None:
		print "Usage: [file to calculate] [name output]"
	else:
		f = open(filename, 'r')
		out = open(outname, 'w')
		filleddict = filldict(f)
		writeoutput(out, filleddict)

if __name__ == '__main__':
	main()
