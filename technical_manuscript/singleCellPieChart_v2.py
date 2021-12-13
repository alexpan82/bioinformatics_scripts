### This script will take any number of clones.txt files generated from the mixcr protocol
### and perform initial analysis on them
### Graphs pie chart of clonal diversity across each sample
### Finds clones shared across all samples (Hopefully becomes multiple venn diagram)
### Tracks clonal frequency across samples (Hopefully)
### Ex: python singleCellPieChart_v2.py [clones.txt files] [output name]

from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
#import matplotlib
#matplotlib.use('Agg')
import pandas as pd
from collections import defaultdict
import numpy
import math

### Write to SAME pdf
#pdf = matplotlib.backends.backend_pdf.PdfPages("%s.pdf" % argv[-1])

#script, paths = argv
# Hold clonal sequences
# dic[clones.txt file name] =[sequence1, freq1, vgene1]
dic = defaultdict(lambda: defaultdict(lambda: []))
uniqseqs = set()
# hold # of clonal seqs
rich = defaultdict(lambda: 0)
# hold count of each clone
countdic = defaultdict(lambda: defaultdict(lambda: 0.0))
#filepaths = open(paths, 'r')
for txtfile in argv[1:len(argv)]:
	txtfile = txtfile.strip()
	key = ["NA"]*8
	value = []
	#print txtfile
	with open(txtfile, 'r') as f:
		for i, line in enumerate(f):
			line = line.strip().split('\t')
			try:
				seq = line[3]
				#freq = float(line[2])
				freq = float(line[1])
				chain = line[5][0:3]
				uniqseqs.add(seq)
				#value.append([seq, freq, chain])
				if i == 1:
					key[0] = seq
					key[1] = line[5].split('*')[0]
					key[2] = line[6].split('*')[0]
					key[3] = line[7].split('*')[0]
					value.append([seq, freq, chain])
					otherfreq =0
				elif i == 2:
					key[4] = seq
					key[5] = line[5].split('*')[0]
					key[6] = line[6].split('*')[0]
					key[7] = line[7].split('*')[0]
					value.append([seq, freq, chain])
					skey = tuple(sorted(key))

				elif i > 2 and i <= 20:
					value.append([seq, freq, chain])
				# Only plot the 1st 20 cdr3 seqs. Thow everything else into 'Others'
				elif i >20:
					otherfreq += freq

			except:
				pass
		if i == 1:
			value.append([seq, 0, chain])
			skey = tuple(sorted(key))
		if i >= 1:
			value.append(['Others', otherfreq, 'IGH'])
			dic[skey][txtfile] = value
			countdic[skey][txtfile] += freq
			del skey
			rich[txtfile] = i
		f.close()

uniqseqs.add('Others')

# Normalize
for s in dic.keys():
	for t in dic[s].keys():
		for n, i in enumerate(dic[s][t]):
			count = dic[s][t][n][1]
			total = countdic[s][t]
			dic[s][t][n][1] = float(count/total)
# Create subplots
#ax = {}

# set font sizes
plt.rcParams.update({'font.size': 18})

# Make a set of sequences (see above) and assign unique color
seqcolors = {}
for seq in uniqseqs:
	seqcolors[seq] = numpy.random.rand(1,3)[0]
seqcolors['Others'] = numpy.array([0.1843, 0.3098, 0.3098])

# Decide how many subplots to make on one pdf page
square = list(map(lambda x: x**2, range(1,21)))
for num in sorted(square):
	if num >= len(argv[1:]):
		size = math.sqrt(num)
		break
#fig = plt.figure()
# Make subplots
#for n in range(0, len(argv[1:])-1):
#	ax[n+1] = fig.add_subplot(1, 1, n+1)

### Pie chart of clonal diversity for each sample
### Need to optimize
i = 0
for key in sorted(dic.keys(), key=lambda x: dic[x].keys(), reverse=False):

 	# Plot Pie charts
        # Normalize colors
	for tfile in dic[key].keys():
		# Plot 'Others' wedge last. So I need to pop it from dic and append on after I've sorted
		fig, ax = plt.subplots()
		print(tfile)
		for t in dic[key][tfile]:
			if t[0] == 'Others':
				othertuple = t
		dic[key][tfile].remove(othertuple)
		# Make sorted frequnecy list (sorted by frequency)
		dvalues = sorted(dic[key][tfile], key=lambda x: x[1], reverse=True)
		sfreqlist = [x[1] for x in dvalues]
		sseqlist = [x[0] for x in dvalues]
		chainlist = [x[2] for x in dvalues]
		# Append 'Others'
		sfreqlist.append(othertuple[1])
		sseqlist.append('Others')
		chainlist.append('IGH')
		# Create patterns list
		patterns = []
		for p in chainlist:
			if p in ['IGK', 'IGL', 'TRB', 'TRD']:
				patterns.append('/')
			else:
				patterns.append('')
		#print '%s\t%s\t%s' % (tfile, key, sfreqlist)
		# Create color list
		scolors=[]
		for s in sseqlist:
			scolors.append(seqcolors[s])
		#print i
		#wedges, texts = ax[i+1].pie(sfreqlist, colors=scolors)
		wedges, texts = ax.pie(sfreqlist, colors=scolors)
		#print sfreqlist
		matplotlib.rcParams['hatch.linewidth'] = 0.4
		matplotlib.rcParams['hatch.color'] = 'k'
		for z, p in enumerate(patterns):
			if p == '/':
				wedges[z].set_hatch(p)
				#wedges[z].set_edgecolor('k')
				#wedges[z].set_linewidth(0.05)

		#ax[i+1].axis('equal')
		ax.axis('equal')
		name = tfile.rstrip('.txt')
		#print name
		#ax[i+1].set_title('%s (%s)' % (name, rich[tfile]))
		ax.set_title('%s (%s)' % (name, rich[tfile]))
		pdf = matplotlib.backends.backend_pdf.PdfPages('%s.pdf' % name)
		#ax[key][i+1].xlabel('Number of unique clones: %s     Frequency of Top Clone: %s' % (len(clnseq[file].values()), max(clnseq[file].values())))
		i += 1

		plt.tight_layout()
		plt.savefig(pdf, format='pdf')
		plt.close()
		pdf.close()

