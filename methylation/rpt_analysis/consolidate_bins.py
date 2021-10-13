# If you have multiple files with this column structure:
# Name \t number
# and some files share some but not all names and you would like 
# to combine these files into one table, look no further!
# This script will do just that

from sys import argv
from collections import defaultdict

dic = defaultdict(lambda: {})

# Write data to dictionary
for i in argv[1:]:
	with open(i,'r') as file:
		for line in file:
			line = line.strip().split('\t')		
			dic[i][line[0]] = line[3]
		file.close()

# Take a union of all the names across the files
names = set()

for i in dic.keys():
	for n in dic[i].keys():
		names.add(n)
	print "\t%s" % i,

for n in sorted(names):
	print "\n%s\t" % n
	for f in dic.keys():
		try:
			print "%s\t" % dic[f][n],
		except:
			print "0\t",

