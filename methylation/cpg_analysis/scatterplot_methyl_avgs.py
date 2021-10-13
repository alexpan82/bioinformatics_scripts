###### Takes output from avg_methyl_reps.py containing cpg information (#chr, pos, cov, numCs, ratio, genesym). First, we identify shared cpgs between all files. Output the avg methylation values of Samples A and B with positional information to graph in a scatterplot. 
from sys import argv
from collections import defaultdict
import sys
script, cpg1, cpg2, genesym, outname = argv

### [pos:[ratio, genesym, cov]]

cpg1_dic = defaultdict( lambda: "-0.01\tNA\tNA" )
cpg2_dic = defaultdict( lambda: "-0.01\tNA\tNA" )

### Hold positional info and find the cpg locations that are shared between both files

#cpg_set = set()
cpg1_set = set()
cpg2_set = set()

output = open(outname, 'w')

with open(cpg1, 'r') as f:
	for line in f:
		line = line.strip()
		if line.split('\t')[5] == genesym:
			#cpg_set.add(line.split('\t')[0] + "\t" + line.split('\t')[1])
			cpg1_set.add(line.split('\t')[0] + "\t" + line.split('\t')[1])
			cpg1_dic[line.split('\t')[0] + "\t" + line.split('\t')[1]] = line.split('\t')[4] + "\t" + line.split('\t')[5] + "\t" + line.split('\t')[2]
	f.close

with open(cpg2, 'r') as f:
	for line in f:
		line = line.strip()
		if line.split('\t')[5] == genesym:
			#cpg_set.add(line.split('\t')[0] + "\t" + line.split('\t')[1])
			cpg2_set.add(line.split('\t')[0] + "\t" + line.split('\t')[1])
			cpg2_dic[line.split('\t')[0] + "\t" + line.split('\t')[1]] = line.split('\t')[4] + "\t" + line.split('\t')[5] + "\t" + line.split('\t')[2] 
	f.close

# Holds positional info for cpgs shared between all files

cpg_total_shared = set.union(cpg1_set, cpg2_set)

output.write("#chr\tpos\tmethyl_value_%s\tmethyl_value_%s\tcov_%s\tcov_%s\tgenesym" % (cpg1, cpg2, cpg1, cpg2))

#for i in cpg_total_shared:
#for i in cpg_set:
for i in cpg_total_shared:
	m_value_cpg1 = float(cpg1_dic[i].split('\t')[0])
	m_value_cpg2 = float(cpg2_dic[i].split('\t')[0])
	if cpg1_dic[i].split('\t')[1] != "NA":
		output.write("\n%s\t%s\t%s\t%s\t%s\t%s" % (i, m_value_cpg1, m_value_cpg2, cpg1_dic[i].split('\t')[2], cpg2_dic[i].split('\t')[2], cpg1_dic[i].split('\t')[1]))
	elif cpg2_dic[i].split('\t')[1] != "NA":
		output.write("\n%s\t%s\t%s\t%s\t%s\t%s" % (i, m_value_cpg1, m_value_cpg2, cpg1_dic[i].split('\t')[2], cpg2_dic[i].split('\t')[2], cpg2_dic[i].split('\t')[1]))











