from sys import argv
from collections import defaultdict

clones = defaultdict(lambda: 0)

with open (argv[1], "r") as f:
	for lines in f:
		linesSplit = lines.strip().split()
		if (linesSplit[2] == "T"):
			

