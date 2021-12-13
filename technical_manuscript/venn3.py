# Expects 4 tables with headers "FUNCTIONAL" and "CDR3" or "CDR3_IMGT"
# Also specify heavy or light chain or both sharing
# Ex: python venn4.py [both|heavy|light] tsv1 tsv2 tsv3 tsv4
 
from sys import argv
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from matplotlib.backends.backend_pdf import PdfPages
import venn
import pandas as pd
from collections import defaultdict


script, chain, f1, f2, f3 = argv

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
		chain_df.append(functional_df[i].loc[~functional_df[i].V_CALL.str.contains('IGH')].copy())

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
# Produce a 4-way venn diagram and output a tsv to std out with this format: 
	# cdr3, # of files with cdr3, files

set_list = []

for i in range(0, len(functional_df)):
	print(i)
	try:
		#tmpdf = set(functional_df[i]['CDR3_IMGT'].tolist())
		#tmpdf = set([x for x in tmpdf if '.' not in x])
		#set_list.append(set(tmpdf))
		set_list.append(set(functional_df[i]['JUNCTION'].tolist()))
	except:
		set_list.append(set(functional_df[i]['CDR3'].tolist()))
		#tmpdf = set(functional_df[i]['CDR3'].tolist())
		#tmpdf = set([x for x in tmpdf if '.' not in x])
		#set_list.append(set(tmpdf))
# Remove '.' from sets
tmp = []
for i in range(0, len(set_list)):
	tmpset = {x for x in set_list[i] if x==x}
	tmpset = set([x for x in tmpset if '.' not in x])
	tmp.append(tmpset)

del set_list
set_list = tmp
del tmp


### Venn diagram of clonal sequences for each sample
pdf = PdfPages('outfile.pdf')
labels = venn.get_labels([set_list[0], set_list[1], set_list[2]],fill=['number','number'])
fig, ax = venn.venn3(labels, names = (argv[2], argv[3], argv[4]))
plt.tight_layout()
plt.savefig(pdf, format='pdf')
plt.close()
pdf.close()


#Now print overlapped sequences
cdr3dic = defaultdict(lambda: [])

for i in range(0, len(set_list)):
	for cdr3 in set_list[i]:
		cdr3dic[cdr3].append(argv[i+2])

for cdr3 in sorted(cdr3dic.keys(), key=lambda x: len(cdr3dic[x]), reverse=True):
	if len(cdr3dic[cdr3]) > 1:
		print('%s\t%s' % (cdr3, '\t'.join(cdr3dic[cdr3])))

