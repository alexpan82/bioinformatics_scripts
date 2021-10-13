# Takes permeth file (#chr\tpos\tcov\tnumCs\tratio) 
# Bins CpGs by cov on a log10 scale

from sys import argv
from collections import defaultdict
import pylab as pl
import numpy as np


script, cpg1 = argv

cpg1_line = open(cpg1, 'r')
cpg1_line.readline()

dic = defaultdict( lambda: 0 )
meth_dic = defaultdict( lambda: [] )
data = []

size = 0
for line in cpg1_line:
	line = line.strip().split('\t')
	cov = int(line[2])
	meth = float(line[3])
	dic[cov] += 1
	meth_dic[cov].append(meth)
	data.append(cov)
	size += 1

x = 10.0
xbins = set()
while x <= 10000:
	x = x * (10.0**0.2)
	xbins.add(x)

initial = 10
print_dic = {}
print_meth = {}
for i in sorted(xbins):
	total = 0
	sum_meth = 0
	len_meth = 0
	totalcov = 0
	for n in dic.keys():
		if initial <= n < i:
			total += dic[n]
			totalcov += n*len(meth_dic[n])
			sum_meth += sum(meth_dic[n])
			#len_meth += len(meth_dic[n])
			
	print_dic[i] = total
	try:
		#print_meth[i] = sum_meth / len_meth
		print_meth[i] = sum_meth / totalcov
	except:
		print_meth[i] = 0
	initial = i
# Uncomment if you want to output bins to display
#
#
print size
print "bin\tnum\tfreq\tavg_meth"
for i in sorted(print_dic.keys()):
	freq = float(print_dic[i])/float(size)
	print "%s\t%s\t%s\t%s" % (i, print_dic[i], freq, print_meth[i])
#
#


# Uncomment if you want to plot bins
#
#
fig, ax = pl.subplots()

xbins.add(10)
xbins = sorted(xbins)
cpg1 = cpg1 + ".pdf"
pl.hist(data, bins=xbins)
pl.gca().set_xscale("log")
pl.xlabel("log10(Coverage per base)", labelpad=30)
pl.ylabel("Number of CpGs")
pl.title("%s" % (cpg1))



ax.annotate("Frequency", xy=(0,-100), xycoords=('figure points', 'data'), va='top', fontsize=10)
ax.annotate("Avg 5mC Ratio", xy=(0,-200), xycoords=('figure points', 'data'), va='top', fontsize=10)

anno_bins = []
print_dic_set = sorted(print_dic.keys())

for i in range(0,len(xbins)-1):
	anno_bins.append((xbins[i] + xbins[i+1])/2)
for s in range(0,len(print_dic_set)):
	freq = float(print_dic[print_dic_set[s]])/float(size)
	freq = round(freq*100.0, 1)
	freq = str(freq)+"%"
	meth = print_meth[print_dic_set[s]]
	meth = round(meth*100.0, 1)
	meth = str(meth)+"%"

	ax.annotate(str(freq), xy=(anno_bins[s],-100), xycoords=('data', 'data'), va='top', ha='center', fontsize=8)
	ax.annotate(str(meth), xy=(anno_bins[s],-200), xycoords=('data', 'data'), va='top', ha='center', fontsize=8)

pl.subplots_adjust(bottom=0.15)

#pl.show()
pl.savefig(cpg1)
#
#
