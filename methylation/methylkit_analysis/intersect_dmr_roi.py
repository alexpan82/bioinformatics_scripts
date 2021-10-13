from sys import argv
import argparse
import csv

# gives the dmc's that overlap with genes and reports them

def main(args):

	# Parse arguments into variables
	dmc = args.dmc
	roi = args.roi
	out = args.out
	dmr = args.dmr
	col = args.col
	if col != 'no':
		col = int(col) - int(1)

	dmc_line = open(dmc, 'r')
	roi_line = open(roi, 'r')
	outfile = open(out, 'w')
	#dmc_line.readline()
	#roi_line.readline()

	dmc_set = set()
	dmc_list = []
	roi_set = set()
	roi_list = []
	out_set = set()
	outfile.write("#chr\tstart\tend\tstrand\tpvalue\tqvalue\tmeth.diff\tgenesym")

	for line in dmc_line:
		dmc_set.add(line.strip())
	for dmc in dmc_set:
	    dmc = dmc.split('\t')
	    dmc_list.append(dmc)


	for line in roi_line:
		roi_set.add(line.strip())
	for roi in roi_set:
		roi = roi.split('\t')
	        roi_list.append(roi)

	if dmr != 'no':
		for rline in roi_list:
			for dline in dmc_list:
				if dline[0] == rline[0]:
					try:
						if int(rline[1])<=int(dline[1])<=int(rline[2]) or int(dline[1])<=int(rline[1])<=int(dline[2]):
							if col != 'no':
								out_set.add("%s\t%s" % ('\t'.join(dline), rline[col]))
								#print "%s\t%s" % ('\t'.join(dline), rline[col])
							else: 
								out_set.add("%s" % '\t'.join(dline))
								#print "%s" % ('\t'.join(dline))
					except:
						pass
	if dmr == 'no':
		for rline in roi_list:
			for dline in dmc_list:
				if dline[0] == rline[0]:
					try:
						if int(rline[1])<=int(dline[1])<=int(rline[2]):
							if col != 'no':
								out_set.add("%s\t%s" % ('\t'.join(dline), rline[col]))
								#print "%s\t%s" % ('\t'.join(dline), rline[col])
							else: 
								out_set.add("%s" % '\t'.join(dline))
								#print "%s" % ('\t'.join(dline))
					except:
						pass
		
	#sort_dmc_set = sorted(out_set)
	for i in out_set:
		outfile.write("\n%s " % i.strip()) 
	dmc_line.close()
	roi_line.close()
	outfile.close()




if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Generates ROI file from a list of genes')
	parser.add_argument('-cpg', dest='dmc', type=str, help='Either DMC/DMR list or .permeth file')
	parser.add_argument('-roi', dest='roi', type=str, help='Region of interest file')
	parser.add_argument('-o', dest='out', type=str, help='Output name')
	parser.add_argument('-r', dest='dmr', default='no', type=str, help='If file inputted was a DMR file or if you want to do a region interesect region analysis, type yes')
	parser.add_argument('-col', dest='col', default='no', type=str, help='If you do not want output with an extra column that contains the genesym associated with each cpg/dmc/dmr input, leave this option alone. If you do, please input a number n that is associated with the nth column in the roi file that contains genesym info')
	args=parser.parse_args()
	main(args)
