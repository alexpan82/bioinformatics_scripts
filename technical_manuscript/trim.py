from sys import argv

h = open (argv[2], "w+")

with open (argv[1], "r") as f:
	i = 0
	for lines in f:
		if i % 2 == 1:
			h.write(lines[:74] + "\n")
		else:
			h.write(lines)
		i += 1

h.close


