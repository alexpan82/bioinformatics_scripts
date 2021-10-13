# Given a tsv and column number (starting with column 1), smooth that column by requiring window_size=1 >= some cutoff. If not, increase window_size and divide new_sum / new window_size
# Ex: python rolling_avg.py -c 2 -t test.tsv
import argparse

def main(args):



	column = int(args.column) - 1
	tsvfile = args.tsv

	# Open tsv and append numbers to list
	data = {'data': []}
	with open(tsvfile, 'r') as f:
		for line in f:
			line = line.strip()
			try:
				line = line.split('\t')[column]
				data['data'].append(float(line))
			except:
				pass
	
	maxbin = max(data['data'])
	# Require cutoff to be > maxbin/100
	cutoff = maxbin/500.0

	# outputlist
	outlist = []
	
	for i, nbin in enumerate(data['data']):
		if nbin >= cutoff:
			outlist.append(nbin)

		else:
			extendback = 0
			extendfwd = 0
			counter = 1
			fwdbins = 0.0
			backbins = 0.0

			while nbin < cutoff:
				# Extending bins may run into edges of list. Need to bound how big the windows are
				if i+counter >= len(data['data'])-1:
					fwd = len(data['data'])-1
				else:
					fwd = i+counter
					fwdbins += 1.0
				if i-counter <= 0:
					back = 0
				else:
					back = i-counter
					backbins += 1.0

				extendback += data['data'][fwd]
				extendfwd += data['data'][back]
				nbin += extendfwd
				nbin += extendback
				counter += 1
				
			# Divide the sum of the bins by the number of bins used
				# counter will 
			outlist.append(nbin/(fwdbins+backbins+1.0))


	for n,i in enumerate(outlist):
		print('%s\t%s' % (n+1,i))





if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='See comments in code')
	parser.add_argument('-c', dest='column', type=str, help='Column in tsv file to take rolling average')
	parser.add_argument('-t', dest='tsv', type=str, help='TSV file')
	args=parser.parse_args()
	main(args)

