# This script makes a 3 bar histogram. If a gene has dmcs that are all in one direction, it gets binned in the (all +) or (all -) bin. If the gene has dmcs in multiple direction, it gets binned in another column
# Takes a dmc/dmr file with meth.diff in column 7 and genesym in column 8
# Only test genes with n number of dmcs

from sys import argv
from collections import defaultdict

script, dmc, n = argv

n = int(n)
# Create dictionary with this structure: {genesym:[]}
dic = defaultdict( lambda: [] )
dmc_set = set()

dmc_line = open(dmc, 'r')

for line in dmc_line:
	dmc_set.add(line.strip())

# If meth.diff is > 0, assign 1 to it. If < 0, assign 0
for dmc in dmc_set:
	dmc = dmc.split('\t')
	try:
		if float(dmc[6]) > 0:
			binary = 1
		elif float(dmc[6]) < 0:
			binary = 0

		dic[dmc[7]].append(binary)
	except:
		pass


hyper = 0
hypo = 0
both = 0

for i in dic.keys():
	if len(dic[i]) >= n:

		if sum(dic[i]) == len(dic[i]):
			hyper += 1
		elif sum(dic[i]) == 0:
			hypo += 1
		elif sum(dic[i]) < len(dic[i]) and sum(dic[i]) > 0:
			both += 1
	else:
		pass

print "hyper\thypo\tboth"
print "%s\t%s\t%s" % (hyper, hypo, both)
