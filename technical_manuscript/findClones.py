from sys import argv
from collections import defaultdict

clones = defaultdict(lambda: 0)
with open (argv[1], "r") as f:
	for lines in f:
		linesStrip = lines.strip().split()
		if (argv[2] in linesStrip[5]):
			clones[linesStrip[1]] += 1

for i in clones.keys():
	print(str(clones[i]) + "\t" + str(i))


