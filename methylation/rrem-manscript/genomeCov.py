# Takes bed file outputted from bedtools genomecov -bg -ibam

from sys import argv

numbases = 0

try:
	covlim = int(argv[2])
except:
	covlim = 0

with open(argv[1], 'r') as f:
	for n, line in enumerate(f):
		line = line.strip().split('\t')
		current_start = int(line[1])
		current_end = int(line[2])
		cov = int(line[3])
		if n == 0:
			prev_end = current_start
		
		if cov >= covlim:
			numbases += current_end - current_start

		if prev_end != current_start:
			numbases += 1

		prev_end = current_end

print(numbases)
