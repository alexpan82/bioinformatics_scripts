from sys import argv

filename = argv[1].replace(".txt", "_removed.txt")

w = open(filename, "w")

with open (argv[1], "r") as f:
	for lines in f:
		if lines.strip().split()[1] == "1.0":
			break
		else:
			w.write(lines)


w.close()
