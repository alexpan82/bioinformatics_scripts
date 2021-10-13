# Summarizes the CDR3 clonal info
# Script takes mixcr clones.txt
# Groups of 3

# Usage: python cdr3BarChart_IGonly_v2.py -m [CLONES.TXT FILES] -c [CHAIN] -n [# ROWS],[# COLS]
# Example: python cdr3BarChart_TRonly_v2.py -m *clones.txt -c heavy -n 4,3

from collections import defaultdict
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter
import argparse
import numpy as np
import math


# Determine whether a given file is mixcr or MakeDb format from the 1st line
def determineFileType(header):
	filetype = ''
	cdr3Index = 0
	chainIndex = 0
	countIndex = 0
	funcIndex = 0
	
	# If 1st line starts with 'cloneID', assume mixcr exportClone output
	if header[0] == 'cloneId':
		filetype = 'mixcr'
		cdr3Index = 3
		chainIndex = 5
		countIndex = 1
		
	elif header[0] == '#count':
		filetype = 'trust'
		cdr3Index = 2
		chainIndex = 4
		countIndex = 0
				
	# Else assume MakeDB output and use JUNCTION as CDR3 sequence
	else:
		try:
			filetype = 'makedb'
			cdr3Index = header.index('JUNCTION')
			chainIndex = header.index('V_CALL')
			funcIndex = header.index('FUNCTIONAL')
			
		except:
			
			print('File type could not be identified. Please input a mixcr, trust, or CHANGEO Makedb formatted file')
			quit()

	return {'filetype': filetype, 'cdr3Index': cdr3Index, 
		'chainIndex': chainIndex, 'countIndex':countIndex, 'funcIndex':funcIndex}
	
# Parse thru line and 
# Read thru mixcr clones.txt or CHANGEO MakeDb tsv and return a dictionary of counts 
def readMixcr(txtlist, uChain):

	# Remove '-'
	filelist = txtlist.copy()
	try:
		filelist.remove('-')
	except:
		pass
	
	# Make sure to only read in lines associated with the correct chain
	chainlist = []
	if uChain == 'heavy':
		chainlist = ['TRB', 'TRD']
	elif uChain == 'light':
		chainlist = ['TRA', 'TRG']
	elif uChain == 'all':
		chainlist = ['TRB', 'TRD', 'TRA', 'TRG']
	else:
		print('Please input either heavy, light, or all for -c')
		quit()
		
	returndic = defaultdict(lambda: defaultdict(lambda: 0))
	cdr3chain = defaultdict(lambda: defaultdict(lambda: ''))
	for txt in filelist:
		with open(txt, 'r') as f:
			
			# Make key name similar to file name
			# name = txt.rstrip('.txt').split('_')[1]
			name = '_'.join(txt.rstrip('.txt').split('_')[0:2])
			
			# Keep track of total reads/cells/contigs in each file
			totalcount = 0.0

			for i, line in enumerate(f):
				line = line.strip().split('\t')
				
				# Determine wheter mixcr or MakeDb format
				if i == 0:
					filetypedic = determineFileType(line)
				
				# Grab CDR3, chaininfo, and counts from proper columns
				#try:
				else:
					try:
						CDR3 = line[filetypedic['cdr3Index']]
						chain = line[filetypedic['chainIndex']][0:3]
						
						if filetypedic['filetype'] == 'makedb':
							count = 1.0
							functional = line[filetypedic['funcIndex']][0]
						else:
							count = float(line[filetypedic['countIndex']])
							functional = 'T'
						
						# Count only if is the IG chain in question
						if (chain in chainlist) and (functional == 'T'):
							returndic[name][CDR3] += count
							totalcount += count
							cdr3chain[name][CDR3] = chain
						
						else:
							pass
	
					except:
						pass

			# Now return percentile values of counts
			for cdr3 in returndic[name].keys():
				returndic[name][cdr3] = returndic[name][cdr3]/totalcount

	return returndic, cdr3chain



def main(args):

	mixcrlist = args.mixcrlist
	cdr3counts, cdr3chain = readMixcr(mixcrlist, args.userChain)

	tmp = []
	for i in mixcrlist:
		try:
			tmp.append('_'.join(i.rstrip('.txt').split('_')[0:2]))
		except:
			tmp.append(i)
	mixcrlist = tmp
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

	# How many suplots to make (do 3 bars in one plot)
	num_rows = int(args.nup.split(',')[0].strip())
	num_cols = int(args.nup.split(',')[1].strip())

	# Plot stacked bars
	# Also make legend
	# legend[gene] = (bar[x][0], frequncy)
	# plt.style.use('ggplot')
	fig, ax = plt.subplots(num_rows, num_cols)

	legendmark = {}

	row_counter = 0
	col_counter = 0
	bin_counter = 0
	bin_size = 3
	for k, f in enumerate(mixcrlist):
		bar = defaultdict(lambda: defaultdict(lambda: None))
		# Color 'Others' bar first if it exists
		if bin_counter == bin_size:
			bin_counter = 0
			col_counter += 1
			if col_counter == num_cols:
				row_counter += 1
				col_counter = 0
		
		
		if 'Others' in plotcdr3counts[f].keys():
			# data = [0] * len(mixcrlist)
			data = [0] * bin_size
			data[bin_counter] = plotcdr3counts[f]['Others']
			bar[f][0] = ax[row_counter, col_counter].bar(np.arange(len(data)), data, width=0.8, linewidth=0.1, color=cdr3colors['Others'])
			legendmark['Others'] = (bar[f][0][0], plotcdr3counts[f]['Others'])
			nextbot = plotcdr3counts[f]['Others']
		else:
			nextbot = 0

		# Color the rest of the stacked bars
		for index, cdr3 in enumerate(sorted(plotcdr3counts[f].keys(), key=lambda u: plotcdr3counts[f][u])):
			# data = [0] * len(mixcrlist)
			data = [0] * bin_size
			data[bin_counter] = plotcdr3counts[f][cdr3]

			if cdr3 != 'Others':
				if cdr3chain[f][cdr3] in ['TRA', 'TRG']:
					if cdr3chain[f][cdr3] == "TRG":
						pattern = 'xxx'
					else:
						pattern = '///'
				else:
					if cdr3chain[f][cdr3] == "TRD":
						pattern = 'ooo'
					else:
						pattern = ''

				bar[f][index+1] = ax[row_counter, col_counter].bar(np.arange(len(data)), data,
						width=0.8, linewidth=0.1, color=cdr3colors[cdr3], bottom=nextbot, hatch=pattern)
				legendmark[cdr3] = (bar[f][index+1][0], plotcdr3counts[f][cdr3])
				nextbot += plotcdr3counts[f][cdr3]

		bin_counter += 1
	# set other bar parameters
	# ax.set_xticks(np.arange(len(mixcrlist)))
	counter = 0
	for i in range(0, num_rows):
		for j in range(0, num_cols):
			ax[i,j].set_xticks(np.arange(bin_size))
			xtickNames = ax[i,j].set_xticklabels(mixcrlist[counter:bin_size+counter])
			plt.setp(xtickNames, rotation=0, fontsize=4)
			ax[i,j].tick_params(axis='y', direction='in', length=1)
			ax[i,j].tick_params(axis='x', length=0)
			ax[i,j].yaxis.set_major_locator(plt.MaxNLocator(5))
			ax[i,j].yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '%s%%' %(int(x*100.0))))
			#ax.set_title('%s (%s%%)     %s (%s%%)' % (name1, mut1, name2, mut2))
			ax[i,j].set_ylim(bottom=0, top=1)
			plt.setp(ax[i,j].get_yticklabels(), fontsize=4)
			counter += bin_size
	

	pdf = PdfPages('bar.pdf')
	plt.tight_layout()
	plt.savefig(pdf, format='pdf')
	plt.savefig('bar.png', dpi=600)
	plt.close()


	pdf.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='See comments in code')
	parser.add_argument('-m', dest='mixcrlist', nargs="+", help="Takes a number of mixcr/changeo/makedb files to be plotted in a series. If you wish to make an empty bar, put '-' in place of file name")
	parser.add_argument('-c', dest='userChain', type=str, help="Please input heavy or light or all")
	parser.add_argument('-n', dest='nup', type=str, help="Shape of subplots: rows,columns")
	args=parser.parse_args()
	main(args)
