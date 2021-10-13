# Takes output from MakeDb.py igblast -i ${output}.fmt7 -s ${input} -r ~osu8725/tools/igblast/database/IMGT/IG[VDJ].fasta --regions --scores --failed
# Finds IGH & IG[KL] pairing info and returns both clonal pairs and summary of pairing

from sys import argv
from collections import defaultdict

# cell_dic[cellname][CDR3] = IG[HKL]
cell_dic = defaultdict(lambda: defaultdict(lambda: 'NA'))

with open(argv[1], 'r') as f:
	for n, line in enumerate(f):
		line = line.strip().split('\t')
		try:
			if n == 0:

				try:
					vcall = line.index('V_CALL')
					func = line.index('FUNCTIONAL')
					try:
						cdr3 = line.index('JUNCTION_10X')
					except:
						cdr3 = line.index('JUNCTION')
					cell = line.index('CELL')
				except:
					print('Table needs to have CELL, V_CALL, FUNCTIONAL, and JUNCTION as column headers')
					quit()

			else:
				if line[func][0] == 'T':
					cellname = line[cell]
					chain = line[vcall][0:3]
					if chain == 'IGK' or chain == 'IGL':
						chain = 'IGLight'
					cdr3seq = chain + ':' + line[cdr3]
					cell_dic[cellname][cdr3seq] = chain
		except:
			pass

handl = 0
horl = 0
other = 0

cdr3_pair = defaultdict(lambda: 0)
hset = set()
lset = set()

for cell in cell_dic.keys():
	chain_list = []
	cdr3_list = []
	for cdr3 in sorted(cell_dic[cell].keys()):
		chain_list.append(cell_dic[cell][cdr3])
		cdr3_list.append(cdr3)
		if cell_dic[cell][cdr3] == 'IGH':
			hset.add(cdr3)
		else:
			lset.add(cdr3)

	#print(chain_list)
	cdr3_pair['\t'.join(cdr3_list)] += 1
	if ('IGH' in chain_list) and ('IGLight' in chain_list) and (len(chain_list) == 2):
		handl += 1
	elif (('IGH' in chain_list) or ('IGLight' in chain_list)) and (len(chain_list) == 1):
		horl += 1
	else:
		other += 1

print('1 HC & 1 LC Pair: %s' % handl)
print('1 HC or 1 LC Only: %s' % horl)
print('Other pairing: %s' % other)
print('Total unique pairings: %s' % len(cdr3_pair.keys()))
print('Total IGH and IGLight CDR3s (respectively): %s\t%s' % (len(hset), len(lset)))

for pair in sorted(cdr3_pair.keys(), key = lambda x: cdr3_pair[x], reverse = True):
	print('%s\t%s' % (cdr3_pair[pair], pair))



	


