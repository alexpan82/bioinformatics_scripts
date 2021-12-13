# Summarizes the mutation / CDR3
# Run after running mixcr for clonal id, trinity to build contigs, and igblast on trinity contigs
# Script takes mixcr clones.txt and igblast/changeo output
# Can take 5 files at a time

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
		filelist = list(filter(('-').__ne__, filelist))
	except:
		pass

	# Make sure to only read in lines associated with the correct chain
	chainlist = []
	if uChain == 'heavy':
		chainlist = ['IGH']
	elif uChain == 'light':
		chainlist = ['IGK', 'IGL']
	elif uChain == 'all':
		chainlist = ['IGK', 'IGL', 'IGH']
	else:
		print('Please input either heavy, light, or all for -c')
		quit()
		
	returndic = defaultdict(lambda: defaultdict(lambda: 0))
	for txt in filelist:
		with open(txt, 'r') as f:
			
			# Make key name similar to file name
			name = txt.rstrip('.txt').split('_')[1]
			
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
						
						else:
							pass
	
					except:
						pass

			# Now return percentile values of counts
			for cdr3 in returndic[name].keys():
				returndic[name][cdr3] = returndic[name][cdr3]/totalcount

	return returndic

# Plot stacked bars
def plotBar(mixcrlist, plotcdr3counts, cdr3colors, axisBreak):
	axisBreak = list(map(lambda x: float(x), axisBreak))
	# 'break' or 'cut-out' the y-axis
	# into two portions - use the top (ax) for the outliers, and the bottom
	# (ax2) for the details of the majority of our data
	fig, (ax, ax2) = plt.subplots(2,1,sharex=True, gridspec_kw={'height_ratios':[1,2]})


	for k, f in enumerate(mixcrlist):
		# Plot same thing in both ax and ax2
		bar = {}
		bar2 = {}
		nextbot = 0
		
		# If empty...
		if len(plotcdr3counts[f].keys()) == 0:
			data = [0] * len(mixcrlist)
			data[k] = 1.
			bar[index+1] = ax.bar(np.arange(len(mixcrlist)), 
					data, width=0.2, linewidth=0.1, 
					color='black', bottom=nextbot)
					
			bar2[index+1] = ax2.bar(np.arange(len(mixcrlist)), 
					data, width=0.2, linewidth=0.1, 
					color='black', bottom=nextbot)			
					
		# Color the rest of the stacked bars
		for index, cdr3 in enumerate(sorted(plotcdr3counts[f].keys(), key=lambda u: plotcdr3counts[f][u])):
			data = [0] * len(mixcrlist)
			data[k] = plotcdr3counts[f][cdr3]

			if cdr3 != 'Others':
				bar[index+1] = ax.bar(np.arange(len(mixcrlist)), 
						data, width=0.7, linewidth=0.1, 
						color=cdr3colors[cdr3], bottom=nextbot)
						
				bar2[index+1] = ax2.bar(np.arange(len(mixcrlist)), 
						data, width=0.7, linewidth=0.1, 
						color=cdr3colors[cdr3], bottom=nextbot)
				nextbot += plotcdr3counts[f][cdr3]
				
		# Color 'Others' bar last if it exists
		if 'Others' in plotcdr3counts[f].keys():
			data = [0] * len(mixcrlist)
			data[k] = plotcdr3counts[f]['Others']
			bar[0] = ax.bar(np.arange(len(data)), data, width=0.7, linewidth=0.1, 
				color='grey', bottom=nextbot)
				
			bar2[0] = ax2.bar(np.arange(len(data)), data, width=0.7, linewidth=0.1, 
				color='grey', bottom=nextbot)




	# Zoom-in / limit the view to different portion of the data
	
	ax.set_ylim(axisBreak[1], 1.)
	ax2.set_ylim(0., axisBreak[0])

	
	# Hide the spines b/t ax and ax2
	ax2.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	#ax.xaxis.tick_bottom()
	#ax2.xaxis.tick_top()
	#ax2.tick_params(labeltop=False)
	
	# set other bar parameters
	ax.set_xticks(np.arange(len(mixcrlist)))
	xtickNames = ax.set_xticklabels(['' if x=='-' else x for x in mixcrlist])
	plt.setp(xtickNames, rotation=0, fontsize=10)
	ax.tick_params(axis='y', direction='in', length=1)
	ax2.tick_params(axis='y', direction='in', length=1)
	ax.tick_params(axis='x', length=0)
	ax2.tick_params(axis='x', length=0, labelsize=24)
	ax2.yaxis.set_major_locator(plt.MaxNLocator(4))
	ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '%s%%' %(int(x*100.0))))
	ax2.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '%s%%' %(int(x*100.0))))
	#ax.set_title('%s (%s%%)     %s (%s%%)' % (name1, mut1, name2, mut2))
	#ax[n+1].axhline(y=0.02)

	pdf = PdfPages('bar.pdf')
	plt.tight_layout()
	plt.savefig(pdf, format='pdf')
	plt.close()
	pdf.close()
	
	

def main(args):

	mixcrlist = args.mixcrlist
	axisBreak = args.axisBreak
	userChain = args.userChain
	cdr3counts = readMixcr(mixcrlist, args.userChain)

	tmp = []
	for i in mixcrlist:
		try:
			tmp.append(i.rstrip('.txt').split('_')[1])
		except:
			tmp.append(i)
	mixcrlist = tmp
	
	allcdr3 = set()
	# Only plot top 10 IGH CDR3's in each file
	plotcdr3counts = defaultdict(lambda: defaultdict(lambda: 0))
	for f in cdr3counts.keys():
		for index, cdr3 in enumerate(sorted(cdr3counts[f].keys(), key=lambda u: cdr3counts[f][u], reverse=True)):
			allcdr3.add(cdr3)

			# This is where the limit of 10 is specified
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
	plotBar(mixcrlist, plotcdr3counts, cdr3colors, axisBreak)
	
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='See comments in code')
	parser.add_argument('-m', dest='mixcrlist', nargs="+", help="Takes a number of mixcr/changeo/makedb files to be plotted in a series. If you wish to make an empty bar, put '-' in place of file name")
	parser.add_argument('-c', dest='userChain', type=str, help="Please input heavy or light or all")
	parser.add_argument('-b', dest='axisBreak', nargs="+", help="Take 2 values to specify where to break y-axis (please input a values b/t 0-1")
	args=parser.parse_args()
	main(args)
