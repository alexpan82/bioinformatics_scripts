###### Takes 3 tab delimenated permeth files (Reps 1, 2, and 3 of Samples A) containing cpg information (#chr, pos, cov, numCs, ratio, genesym). Then figure out which cpgs are shared between all files. If a cpg has above 10 read cov in any of the 3 files, retain that cpg for all files. Then, we take the avg methylation value (ratio) = total Cs/total cov.
###### The above functionallity has been commented out; now this script only retains those cpgs that are >= 10 cov in all 3 files
from collections import defaultdict
from sys import argv
script, cpg1, cpg2, cpg3, outname = argv


### [pos:[cov, numCs, genesym]]
cpg1_dic = defaultdict( lambda: "0\t0\tNA" )
cpg2_dic = defaultdict( lambda: "0\t0\tNA" )
cpg3_dic = defaultdict( lambda: "0\t0\tNA" )

### Hold positional info and find the cpg locations that are shared between all 3 files
cpg1_set = set()
cpg2_set = set()
cpg3_set = set()

output = open(outname, 'w')

with open(cpg1, 'r') as f:
	for line in f:
		line = line.strip().split('\t')
		cpg1_set.add(line[0] + "\t" + line[1])
		try:
			cpg1_dic[line[0] + "\t" + line[1]] = line[2] + "\t" + line[3] + "\t" + line[5]
		except:
			cpg1_dic[line[0] + "\t" + line[1]] = line[2] + "\t" + line[3] + "\t" + "hi"
	f.close

with open(cpg2, 'r') as f:
	for line in f:
		line = line.strip().split('\t')
		cpg2_set.add(line[0] + "\t" + line[1])
		try:

			cpg2_dic[line[0] + "\t" + line[1]] = line[2] + "\t" + line[3] + "\t" + line[5]
		except:
			cpg2_dic[line[0] + "\t" + line[1]] = line[2] + "\t" + line[3] + "\t" + "hi"
	f.close

with open(cpg3, 'r') as f:
	for line in f:
		line = line.strip().split('\t')
		cpg3_set.add(line[0] + "\t" + line[1])
		try:
			cpg3_dic[line[0] + "\t" + line[1]] = line[2] + "\t" + line[3] + "\t" + line[5]
		except:
			cpg3_dic[line[0] + "\t" + line[1]] = line[2] + "\t" + line[3] + "\t" + "hi"
	f.close

# Holds positional info for cpgs shared between all files
cpg_total_shared = set.union(cpg1_set, cpg2_set, cpg3_set)

print "%s\t%s\t%s\t%s" % (len(cpg1_set), len(cpg2_set), len(cpg3_set), len(cpg_total_shared))
### Apply 10 read cov filter. Again, if one file has > 10, then retain that cpg across all files

output.write("#chr\tpos\ttotal_cov\ttotal_numCs\tratio\tgenesym")

count = 0

for i in cpg_total_shared:
	try:
		temp1 = cpg1_dic[i].split('\t')
		temp2 = cpg2_dic[i].split('\t')
		temp3 = cpg3_dic[i].split('\t')
		c1 = int(temp1[0])
		c2 = int(temp2[0])
		c3 = int(temp3[0])
		nc1 = int(temp1[1])
		nc2 = int(temp2[1])
		nc3 = int(temp3[1])		

		#if int(cpg1_dic[i].split('\t')[0]) >= 10 or int(cpg2_dic[i].split('\t')[0]) >= 10 or int(cpg3_dic[i].split('\t')[0]) >= 10:
		if c1 >= 10 and c2 >= 10 and c3 >= 10:
			count += 1
			total_cov = c1 + c2 + c3
			total_numc = nc1 + nc2 + nc3
			total_ratio = float(total_numc) / float(total_cov)
			if temp1[2] != "NA":
				if temp1[2] == "hi":
					output.write("\n%s\t%s\t%s\t%s" % (i, total_cov, total_numc, total_ratio))
				else:
					output.write("\n%s\t%s\t%s\t%s\t%s" % (i, total_cov, total_numc, total_ratio, temp1[2]))
			elif temp2[2] != "NA":
				if temp2[2] == "hi":
					output.write("\n%s\t%s\t%s\t%s" % (i, total_cov, total_numc, total_ratio))
				else:
					output.write("\n%s\t%s\t%s\t%s\t%s" % (i, total_cov, total_numc, total_ratio, temp2[2]))
			elif temp3[2] != "NA":
				if temp3[2] == "hi":
					output.write("\n%s\t%s\t%s\t%s" % (i, total_cov, total_numc, total_ratio))
				else:
					output.write("\n%s\t%s\t%s\t%s\t%s" % (i, total_cov, total_numc, total_ratio, temp3[2]))

	except:
		pass
print count
