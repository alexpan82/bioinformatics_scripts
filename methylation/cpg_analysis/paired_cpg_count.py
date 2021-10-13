# This script identifies CpGs that are from the same CpG dinucleotide pair (top and bottom strand CpGs). Then it outputs statistics such as mean methylation difference between top and bottom CpGs, std devs of top and bottom, etc.  

# Outputs 2 histograms. One shows the number of CpG pairs vs their methylation difference
# The second shows number of CpG pairs vs their coverage difference

# For now, just takes a permeth file. Outputs the same permeth file but with * in the last column denoting pairs

from sys import argv
from collections import defaultdict
import numpy as np
from scipy import stats
import csv
import matplotlib.pyplot as pl
import seaborn as sns
import math

script, permeth = argv

# Declare a dictionary. Has this structure:
# {chr:{position:(coverage, numC, ratio)}}

dic = defaultdict( lambda: {} )

# Read from file and set up dictionary format
with open(permeth) as tsvfile:
	tsvreader = csv.reader(tsvfile, delimiter="\t")
	for line in tsvreader:
		if line[0].startswith("c"):
			dic[line[0]][int(line[1])] = (int(line[2]), int(line[3]), float(line[4]))

	tsvfile.close()

# Find paired CpGs and put into a dictionary of sets
# To do this, the script checks whether the current position - 1 = the last position
# If so, then we know the current and last position are "paired"
# Then, the last position = -1 so that the next current position - 1 will never equal the last position
# If current pos != last pos, then the last_pos will = current pos at the end of the loop
# {chr:{position:(coverage, numC, ratio)}}
last_pos = -1
diffmeth_list = []
diffcov_list = []
top_meth = ()
bot_meth = ()
top_cov = ()
bot_cov = ()
top_pos = ()
bot_pos = ()

count = 0

for chrm in dic.keys():
	for pos in sorted(dic[chrm].keys()):
		if pos - 1 == last_pos:
			top_meth = top_meth + (dic[chrm][last_pos][2],)
			bot_meth = bot_meth + (dic[chrm][pos][2],)
			top_cov = top_cov + (dic[chrm][last_pos][0],)
			bot_cov = bot_cov + (dic[chrm][pos][0],)
			top_pos = top_pos + (last_pos,)
			bot_pos = bot_pos + (pos,)
			# Top - Bottom
			diffmeth = dic[chrm][last_pos][2] - dic[chrm][pos][2]
			diffmeth_list.append(diffmeth)
			diffcov = dic[chrm][last_pos][0] - dic[chrm][pos][0]
			diffcov_list.append(diffcov)
			if diffmeth > 0.3 and diffmeth < 0.7:
				count += 1
			# add a demarcation to the original dictionary so when we rewrite, we will
			# be able to tell which lines contain a CpG pair
			dic[chrm][last_pos] = dic[chrm][last_pos] + ("*",)
			dic[chrm][pos] = dic[chrm][pos] + ("*",)
			last_pos = -1
		else:
			last_pos = pos 


abs_diffmeth_list = map(abs, diffmeth_list)
abs_diffcov_list = map(abs, diffcov_list)

print("Number of CpG pairs:\t%s" % len(diffmeth_list))
print("Average of the absolute 5mC differences between CpG pairs:\t%s" % np.mean(abs_diffmeth_list))
print("Average of the absolute coverage differences between CpG pairs:\t%s" % np.mean(abs_diffcov_list))
print("Standard deviation of the absolute 5mC differences between CpG pairs:\t%s" % np.std(abs_diffmeth_list))
print("Standard deviation of the absolute coverage differences between CpG pairs:\t%s" % np.std(abs_diffcov_list))
print("Paired T-test P-value between Top and Bottom Strand CpG Methylation:\t%s" % stats.ttest_rel(top_meth, bot_meth)[1])
print("Paired T-test P-value between Top and Bottom Strand coverage:\t%s" % stats.ttest_rel(top_cov, bot_cov)[1])

print("Number of CpG pairs that have between 0.3 and 0.7 methylation difference:\t%s" % count)


outfile = open(permeth, 'w')

for chrm in dic.keys():
	for pos in sorted(dic[chrm].keys()):
		outfile.write("%s\t%s\t%s\n" % (chrm, pos, '\t'.join(map(str, dic[chrm][pos]))))

outfile.close()


# Make pdf output names and define bins of histogram
outname1 = permeth + ".meth_diff.pdf"
outname2 = permeth + ".cov.pdf"
#outname3 = permeth + ""


# Plot density plots
sns.set()
plot = sns.kdeplot(diffmeth_list, bw=0.01, shade=True)
pl.xlabel("Methylation difference (Top - Bottom)")
pl.ylabel("Number of CpGs Dinucleotide Pairs")
pl.title("%s" % (outname1))
pl.savefig(outname1)

pl.close()
plot2 = sns.kdeplot(diffcov_list, bw=10, shade=True)
pl.xlabel("Coverage difference (Top - Bottom)")
pl.ylabel("Number of CpGs Dinucleotide Pairs")
pl.title("%s" % (outname2))
pl.savefig(outname2)
pl.close()

# Plot line
pl.plot(top_pos, top_meth, color='blue', label='Top CpG Methylation')
pl.plot(bot_pos, bot_meth, color='red', label='Bottom CpG Methylation')
pl.xlabel("Position")
pl.ylabel("Methylation")
#pl.xticks(fontsize=6)
ax = pl.gca()
ax.ticklabel_format(useOffset=False, style='plain')
pl.title("CpG Dinucleotide Pair Methylation")
pl.legend()
pl.show()


