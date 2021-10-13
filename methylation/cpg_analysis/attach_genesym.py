from sys import argv
# The 4 arguments are a list with #chr \t position, a bed file (/fs/scratch/ccri0063/amal/roi/generate_roi.py), a region name (just input string), and out name

# attaches genesym to lines in 1st arg

script, dmc, roi, region, out = argv

dmc_line = open(dmc, 'r')
roi_line = open(roi, 'r')
header = dmc_line.readline().strip()
header = "%s\t%s" % (header, region)
num = len(header.split('\t'))
roi_line.readline()
dmc_dic = {}
dmc_set = set()
roi_dic = {}
roi_set = set()
dmc_set = set()

# Skip header
iterdmc = iter(dmc_line)
for line in iterdmc:
	#line = line.strip()
	linex = line.strip().split('\t')
	chrm = linex[0]
	position = int(linex[1])
	if len(linex) <= num:
		tab = []
		diff = num - len(linex) - 1
		#print diff
		for i in range(0, diff):
			tab.append('\t')
	#print '\t'.join(linex)
	try:
		dmc_dic[chrm][position] = '\t'.join(linex) + ''.join(tab)
	except:
		dmc_dic[chrm] = {}
		dmc_dic[chrm][position] = '\t'.join(linex) + ''.join(tab)

iterroi = iter(roi_line)
for line in iterroi:
	linex = line.strip().split('\t')
	try:
		roi_set.add("%s\t%s\t%s\t%s" % (linex[0],linex[1],linex[2],linex[3]))
	except:
		roi_set.add("%s\t%s\t%s\tintergenic" % (linex[0],linex[1],linex[2]))
for line in roi_set:
	line = line.strip().split('\t')
	chrm = line[0]
	rpos = (int(line[1]), int(line[2]))
	try:
		roi_dic[chrm][rpos] = line[3]
	except:
		roi_dic[chrm] = {}
		roi_dic[chrm][rpos] = line[3]
dmc_line.close()
roi_line.close()

outfile = open(out, 'w')
outfile.write(header)
num = len(header.split('\t'))
assDmc = {}
#test = set()
#test.add("chr7")
for c in dmc_dic.keys():#test:#
	#print c
	assDmc[c] = {}
	for dpos in dmc_dic[c].keys():
		for rrange in roi_dic[c]:
			if dpos >= rrange[0] and dpos <= rrange[1]:
				try:
					#print dmc_dic[c][dpos]
					assDmc[c][dpos][roi_dic[c][rrange]] = "%s\t%s" % (dmc_dic[c][dpos], roi_dic[c][rrange])
				except:
					assDmc[c][dpos] = {}
					assDmc[c][dpos][roi_dic[c][rrange]] = "%s\t%s" % (dmc_dic[c][dpos], roi_dic[c][rrange])
#sort_dmc_set = sorted(dmc_set)
notAssDmc = {}
count = 0
rcount = 0
dcount = 0
for c in assDmc.keys():
	notAssDmc[c] = set(dmc_dic[c].keys()) - set(assDmc[c].keys())
	count += len(set(dmc_dic[c].keys()))
	dcount += len(set(assDmc[c].keys()))
	#print set(assDmc[c].keys())
	#print len(notAssDmc[c])
	for pos in assDmc[c].keys():
		rcount += len(assDmc[c][pos].keys())
		for r in assDmc[c][pos].keys():
			outfile.write("\n%s" % (assDmc[c][pos][r]))
	for pos in notAssDmc[c]:
		outfile.write("\n%s" % dmc_dic[c][pos])
print "Given %s and %s" % (dmc, region)
print "Total number of unique DMCs: %s" % count
print "Number of regions with DMCs: %s" % rcount
print "Number of unique DMCs in all regions: %s" % dcount
