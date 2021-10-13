# Summarizes the mutation / CDR3
# Run after running mixcr for clonal id, trinity to build contigs, and igblast on trinity contigs
# Script takes mixcr clones.txt and igblast/changeo output
# Only takes 3 timepoints at a time

# Usage: python mutation_Bars.py -m A_clones.txt B_clones.txt C_clones.txt -i A_igblast.tsv B_igblast.tsv C_igblast.tsv
# Usage: python mutation_Bars.py -m A_clones.txt B_clones.txt - -i A_igblast.tsv B_igblast.tsv -

from collections import defaultdict
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter
import argparse
import numpy as np

# Read thru mixcr clones.txt file and return a dictionary of counts 
def readMixcr(txtlist):
	returndic = defaultdict(lambda: defaultdict(lambda: 0))
	for txt in txtlist:
		try:
			with open(txt, 'r') as f:

				totalcount = 0

				for line in f:
					line = line.strip().split('\t')
					try:
						CDR3 = line[3]
						chain = line[5][0:3]
						count = float(line[2])
						if chain == 'IGH':
							returndic[txt][CDR3] += count
							totalcount += count
					except:
						pass

				# Now return percentile values of counts
				for cdr3 in returndic[txt].keys():
					returndic[txt][cdr3] = returndic[txt][cdr3]/totalcount
		except:
			pass
	return returndic

# Read thru igblast file and return a dictionary of mutation rates and vgene lengths
def readIG(txtlist):
	returndic = defaultdict(lambda: defaultdict(lambda: [[0,0]]))

	for txt in txtlist:
		try:
			with open(txt, 'r') as f:
				for line in f:
					line = line.strip().split('\t')
					try:
						CDR3 = line[28]
						mutation = 1.0 - float(line[32])
						vlen = int(line[13])
						functional = line[2]
						chain = line[7][0:3]

						if chain == 'IGH' and functional == 'T':
							returndic[txt][CDR3].append([mutation, vlen])
							returndic[txt][CDR3].sort(key=lambda x: x[1], reverse=True)
					except:
						pass
		except:
			pass


	return returndic

def main(args):

	cdr3counts = readMixcr(args.mixcrlist)
	mutationdic = readIG(args.iglist)

	allcdr3 = set()
	# Only plot top 10 IGH CDR3's in each file
	plotcdr3counts = defaultdict(lambda: defaultdict(lambda: 0))
	for f in cdr3counts.keys():
		for index, cdr3 in enumerate(sorted(cdr3counts[f].keys(), key=lambda u: cdr3counts[f][u], reverse=True)):
			allcdr3.add(cdr3)

			if index > 9:
				plotcdr3counts[f]['Others'] += cdr3counts[f][cdr3]
			else:
				plotcdr3counts[f][cdr3] = cdr3counts[f][cdr3]


	# Assign color to each cdr3
	cdr3colors = defaultdict(lambda: [0,0,0])
	for cdr3 in allcdr3:
		cdr3colors[cdr3] = np.random.rand(1,3)[0]
	# Make 'Others' gray
	cdr3colors['Others'] = np.array([0.1843, 0.3098, 0.3098])


	# Plot stacked bars
	# Also make legend
	# legend[gene] = (bar[x][0], frequncy)
	fig, ax = plt.subplots()

	legendmark = {}

	for k, f in enumerate(args.mixcrlist):
		bar = {}
		# Color 'Others' bar first if it exists
		if 'Others' in plotcdr3counts[f].keys():
			data = [0] * len(args.mixcrlist)
			data[k] = plotcdr3counts[f]['Others']
			bar[0] = ax.bar(np.arange(len(data)), data, width=0.5, linewidth=0.1, color=cdr3colors['Others'])
			legendmark['Others'] = (bar[0][0], plotcdr3counts[f]['Others'])
			nextbot = plotcdr3counts[f]['Others']
		else:
			nextbot = 0

		# Color the rest of the stacked bars
		for index, cdr3 in enumerate(sorted(plotcdr3counts[f].keys(), key=lambda u: plotcdr3counts[f][u])):
			data = [0] * len(args.mixcrlist)
			data[k] = plotcdr3counts[f][cdr3]

			if cdr3 != 'Others':
				bar[index+1] = ax.bar(np.arange(len(args.mixcrlist)), data, width=0.5, linewidth=0.1, color=cdr3colors[cdr3], bottom=nextbot)
				legendmark[cdr3] = (bar[index+1][0], plotcdr3counts[f][cdr3])
				nextbot += plotcdr3counts[f][cdr3]


	# set other bar parameters
	ax.set_xticks(np.arange(len(args.mixcrlist)))
	xtickNames = ax.set_xticklabels(args.mixcrlist)
	plt.setp(xtickNames, rotation=0, fontsize=4)
	ax.tick_params(axis='y', direction='in', length=1)
	ax.tick_params(axis='x', length=0)
	ax.yaxis.set_major_locator(plt.MaxNLocator(5))
	ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '%s%%' %(int(x*100.0))))
	#ax.set_title('%s (%s%%)     %s (%s%%)' % (name1, mut1, name2, mut2))
	ax.set_ylim(bottom=0, top=1)
	#ax[n+1].axhline(y=0.02)


	# Make legend
	legendtuple = tuple()
	nametuple = tuple()
	for cdr3 in sorted(legendmark.keys(), key=lambda u: legendmark[u][1], reverse=True):
		if cdr3 != 'Others':
			legendtuple += (legendmark[cdr3][0],)
			# add mutation and vlen info to legend
			# only add the mutation that has the longest vlen
			nameinfo = ''
			for i in args.iglist:
				vlen = mutationdic[i][cdr3][0][1]
				mutation = mutationdic[i][cdr3][0][0]
				addon = ''

				if vlen < 290 and vlen > 200:
					addon = '*'
				elif vlen < 200:
					addon = '**'

				if vlen > 0:
					tmp = str(round(mutation*100., 2)) + '%' + addon
				elif vlen == 0:
					tmp = 'NA'

				nameinfo += tmp
				nameinfo += ', '
			nametuple += (nameinfo,)

	if 'Others' in legendmark.keys():
		legendtuple += (legendmark['Others'][0],)
		nametuple += ('Others',)
	ax.legend(legendtuple, nametuple, loc=0, fontsize=8)

	

	pdf = PdfPages('bar.pdf')
	plt.tight_layout()
	plt.savefig(pdf, format='pdf')
	plt.close()


	pdf.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='See comments in code')
	parser.add_argument('-m', dest='mixcrlist', nargs="+", help="Takes 3 mixcr clones.txt files to be plotted in a series. If < 3 files, put '-' in place of file name")
	parser.add_argument('-i', dest='iglist', nargs="+", help="Takes 3 igblast/changeo tsv files to be plotted in a series. Must correspond the order that mixcr files are given. If < 3 files, put '-' in place of file name")
	args=parser.parse_args()
	main(args)