###### % methyation for every methylated CpG with read cov >= 10 in a given .roi for a .permeth
###### 1st arg is .permeth file, 2nd is an roi (cols are #chr \t start \t end \t genename), 3rd is short output identifier

from sys import argv
script, permeth, roi, out = argv

permeth_line = open(permeth, 'r')
roi_line = open(roi, 'r')

autophagy_set = set()
roi_line_list = []
for line in roi_line:
	line = line.strip()
	autophagy_set.add(line.split('\t')[3])
	roi_line_list.append(line)

for pline in permeth_line:
	pline = pline.strip().split('\t')
	try:
		if int(pline[2]) >= 10:
			for rline in roi_line_list:
				rline = rline.strip().split('\t')
				if pline[0] == rline[0]:
					try:
						if int(rline[1])<=int(pline[1])<=int(rline[2]):
							try:
								with open(out, 'a') as f:
									f.write("%s\t%s\n" % (pline, rline[3]))
							except:
								with open(out, 'a') as f:
									f.write("%s\n" % pline)
					except:
						pass
	except:
		pass
permeth_line.close()
roi_line.close()
