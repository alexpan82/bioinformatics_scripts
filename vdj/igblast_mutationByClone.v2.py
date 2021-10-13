# Takes output from MakeDb.py igblast -i ${output}.fmt7 -s ${input} -r ~osu8725/tools/igblast/database/IMGT/IG[VDJ].fasta --regions --scores --failed
# Calculates mutation rate by CDR3 clone

from sys import argv
from collections import defaultdict

# vmutation[cdr3] = [total matches, total v length, vgenename, totalcontigs, set(cell barcodes), set(constant region)]
vmutation = defaultdict(lambda: [0,0,defaultdict(lambda: 0),0, set(), defaultdict(lambda: 0)])

with open(argv[1], 'r') as f:
	for n, line in enumerate(f):
		line = line.strip().split('\t')
		try:
			if n == 0:
				vcall = line.index('V_CALL')
				func = line.index('FUNCTIONAL')
				seqvlength = line.index('V_SEQ_LENGTH')
				v_id = line.index('V_IDENTITY')
				cdr3 = line.index('JUNCTION')
				#cdr3 = line.index('CDR3_IMGT')
				try:
					cell = line.index('CELL')
					c_call = line.index('C_CALL')
				except:
					cell = line.index('SEQUENCE_ID')
					c_call = 1000
				#vmatches = round(v_id * seqvlength, 0)
				continue
			if line[func] == 'T':
				cdr3seq = line[cdr3]
				vmatches = round(float(line[v_id]) * float(line[seqvlength]), 0)
				vmutation[cdr3seq][0] += vmatches
				vmutation[cdr3seq][1] += float(line[seqvlength])
				#vmutation[cdr3seq][2] = line[vcall] + ' ' + line[vcall+1] + ' ' + line[vcall+2]
				vmutation[cdr3seq][2][line[vcall] + ' ' + line[vcall+1] + ' ' + line[vcall+2]] += 1
				vmutation[cdr3seq][3] += 1
				vmutation[cdr3seq][4].add(line[cell])
				if line[c_call] == '' or c_call == 1000:
					vmutation[cdr3seq][5]['NA'] += 1
				else:
					vmutation[cdr3seq][5][line[c_call]] += 1
		except:
			pass
			#sampleid = line[0].split('_')[0]
			#func = line[2]
			#seqvlength = line[13]
			#v_id = line[32]
			#cdr3 = line[-2]
			#print(sampleid, func, seqvlength, cdr3, v_id)

if len(argv) == 3:
	specifiedCDR3 = set()
	with open(argv[2], 'r') as f:
		for line in f:
			line = line.strip()
			specifiedCDR3.add(line)


	for cdr3 in sorted(specifiedCDR3):
		print('%s\t%s\t%s\t%s\t%s\t%s' % (cdr3, vmutation[cdr3][2], (1.0-vmutation[cdr3][0]/vmutation[cdr3][1]), vmutation[cdr3][3], len(vmutation[cdr3][4]), ', '.join(list(vmutation[cdr3][5]))))

else:


	for cdr3 in sorted(vmutation.keys(), key=lambda x: vmutation[x][3], reverse=True):
		c_call_list = []
		for c in sorted(vmutation[cdr3][5], key=lambda x: vmutation[cdr3][5][x], reverse=True):
			c_call_list.append('%s (%s)' % (c, vmutation[cdr3][5][c]))
		align_list = []
		for a in sorted(vmutation[cdr3][2], key=lambda x: vmutation[cdr3][2][x], reverse=True):
			align_list.append('%s (%s)' % (a, vmutation[cdr3][2][a]))
		print('%s\t%s\t%s\t%s\t%s\t%s' % (cdr3, '; '.join(align_list), (1.0-vmutation[cdr3][0]/vmutation[cdr3][1]), vmutation[cdr3][3], len(vmutation[cdr3][4]), ', '.join(c_call_list)))
