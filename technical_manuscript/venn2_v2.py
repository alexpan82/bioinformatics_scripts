# Expects 4 tables with headers "FUNCTIONAL" and "CDR3" or "CDR3_IMGT"
# Also specify heavy or light chain or both sharing
# Ex: python venn4.py [both|heavy|light] tsv1 tsv2 tsv3 tsv4
 
from sys import argv
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
import pandas as pd
from collections import defaultdict


script, chain, f1, f2 = argv

# Write to data frame
df = []
for f in argv[2:]:
	df.append(pd.read_csv(f, delimiter='\t', index_col=0))

# Parse mixcr output: replace targetSequences and allVHitsWithScore with JUNCTION and V_CALL
for i in range(0, len(df)):
	if "targetSequences" in df[i].columns:
		df[i] = df[i].rename(columns={"targetSequences":"JUNCTION", "allVHitsWithScore":"V_CALL"})

functional_df = []

# Filter by functional cdr3's only
for i in range(0, len(df)):
	try:
		functional_df.append(df[i].loc[df[i].FUNCTIONAL == 'T'].copy())
	except:
		try:
			functional_df.append(df[i].loc[df[i].FUNCTIONAL == True].copy())
		except:
			print('%s does not have a column header named FUNCTIONAL. Will assume everything is functional' % i)
			functional_df.append(df[i])

del df

# Filter by chain
if argv[1] == 'light':
	chain_df = []
	for i in range(0, len(functional_df)):
		#chain_df.append(functional_df[i].loc[~functional_df[i].V_CALL.str.contains('IGH')].copy())
		chain_df.append(functional_df[i].loc[functional_df[i].V_CALL.str.contains('IGK') | functional_df[i].V_CALL.str.contains('IGL')].copy())

	del functional_df
	functional_df = chain_df
	del chain_df

elif argv[1] == 'heavy':
	chain_df = []
	for i in range(0, len(functional_df)):
		chain_df.append(functional_df[i].loc[functional_df[i].V_CALL.str.contains('IGH')].copy())

	del functional_df
	functional_df = chain_df
	del chain_df

else:
	pass

# Pass CDR3 to sets and intersect

set_list = []

for i in range(0, len(functional_df)):
	try:
		set_list.append(set(functional_df[i]['JUNCTION'].tolist()))
	except:
		set_list.append(set(functional_df[i]['CDR3'].tolist()))

# Remove '.' from sets
tmp = []
for i in range(0, len(set_list)):
	tmpset = {x for x in set_list[i] if x==x}
	tmpset = set([x for x in tmpset if '.' or '' not in x])
	tmp.append(tmpset)

del set_list
set_list = tmp
del tmp

# Overlap
sharedCDR3 = set_list[0] & set_list[1]
intersect = len(sharedCDR3)

# Add up # of reads/contigs corresponding to each chain
totals = [0,0]

for i in range(0, len(functional_df)):
	try:
		totals[i] += functional_df[i]['cloneCount'].sum()
	except:
		totals[i] += len(functional_df[i])

# Add up # of reads/contigs corresponding to shared clones
sharedcount = [0,0]
for i in range(0, len(functional_df)):
	try:
		sharedcount[i] += functional_df[i].loc[functional_df[i].JUNCTION.isin(sharedCDR3)]['cloneCount'].sum()
	except:
		sharedcount[i] += len(functional_df[i].loc[functional_df[i].JUNCTION.isin(sharedCDR3)])

print(totals, sharedcount)

### Venn diagram of clonal sequences for each sample
fig = plt.figure()
ax = fig.add_subplot(111)

#colors = {'10X': 'skyblue', 'Immunoseq': 'purple', 'lcRNA-Seq': 'red', 'mRNA-Seq': 'green', '1ng': 'orange'}
colors = defaultdict(lambda: 'white')

abundance1 = int(100.*sharedcount[0]/totals[0])
abundance2 = int(100.*sharedcount[1]/totals[1])
lowest_abun = min([abundance1, abundance2])
intersect_label = "(%s%%/%s%%)" % (abundance1, abundance2)
print('hi')
print((abundance1-(lowest_abun/2), abundance2-(lowest_abun/2), lowest_abun/2))
figure = venn2(subsets = (abundance1-(lowest_abun/2), abundance2-(lowest_abun/2), lowest_abun/2), 
	set_labels = (argv[2].split('.')[0], argv[3].split('.')[0]),
	set_colors=(colors[argv[2].split('.')[0]], colors[argv[3].split('.')[0]]))
venn2_circles(subsets = (abundance1-(lowest_abun/2), abundance2-(lowest_abun/2), lowest_abun/2),
	linewidth=2)
print(len(set_list[0])-intersect, len(set_list[1])-intersect, intersect)

for text in figure.set_labels:
	text.set_fontsize(36)
for text in figure.subset_labels:
	text.set_fontsize(36)

figure.get_label_by_id('10').set_text('%s' % str(len(set_list[0])-intersect))
figure.get_label_by_id('01').set_text('%s' % str(len(set_list[1])-intersect))
figure.get_label_by_id('11').set_text('%s' % str(intersect))

lbl1 = figure.get_label_by_id("A")
lbl2 = figure.get_label_by_id("B")

x1, y1 = lbl1.get_position()
x2, y2 = lbl2.get_position()

lbl1.set_position((x1-0.02, y1+0.25))
lbl2.set_position((x2+0.02, y2+0.25))

print(x1,y1,x2,y2)

ax.text(0.5, 0.5-0.08, intersect_label, transform=ax.transAxes, fontsize=28, verticalalignment='top', horizontalalignment='center')
# move label


'''
def label_by_id(label, ID):
	if label.endswith('\n') is False:
		label += '\n'
	tmplabel = 5*['']
	num = figure.get_label_by_id(ID).get_text()
	tmplabel[2] = num
	tmplabel[4] = label

	label = '\n'.join(tmplabel)
	figure.get_label_by_id(ID).set_text(label)

intersect_label = "(%s%%/\n%s%%)" % (int(100.*sharedcount[0]/totals[0]), int(100.*sharedcount[1]/totals[1]))
labels = [argv[2].split('.')[0], argv[3].split('.')[0], intersect_label]
for label, ID in zip(labels, ['10','01','11']):
	label_by_id(label, ID)
'''
outname = argv[2].split('.')[0] + "_" + argv[3].split('.')[0] + "_" + argv[1] + ".png"
plt.tight_layout()
plt.savefig(outname, format='png')
plt.close()
#pdf.close()


#Now print overlapped sequences
"""
cdr3dic = defaultdict(lambda: [])

for i in range(0, len(set_list)):
	for cdr3 in set_list[i]:
		cdr3dic[cdr3].append(argv[i+2])

for cdr3 in sorted(cdr3dic.keys(), key=lambda x: len(cdr3dic[x]), reverse=True):
	if len(cdr3dic[cdr3]) > 1:
		print('%s\t%s' % (cdr3, '\t'.join(cdr3dic[cdr3])))

"""
