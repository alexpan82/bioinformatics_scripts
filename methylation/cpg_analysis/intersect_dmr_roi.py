from sys import argv
# The 3 arguments are the dmc list from methylkit's calculateDiffMeth(), roi file (/fs/scratch/ccri0063/amal/roi/generate_roi.py), and output file name
# gives the dmc's that overlap with genes and reports them

script, dmc, roi, out = argv

dmc_line = open(dmc, 'r')
roi_line = open(roi, 'r')
outfile = open(out, 'w')
dmc_line.readline()
roi_line.readline()

dmc_set = set()
roi_set = set()
out_set = set()
outfile.write("#chr\tpos\tcov\tnumC\t5mCratio\tgenesym")

for line in dmc_line:
	dmc_set.add(line.strip())
for line in roi_line:
	roi_set.add(line.strip())

for rline in roi_set:
	rlinex = rline.strip().split('\t')
	for dline in dmc_set:
		dlinex = dline.strip().split('\t')
		if dlinex[0] == rlinex[0]:
			if int(rlinex[1])<=int(dlinex[1])<=int(rlinex[2]) or int(dlinex[1])<=int(rlinex[1])<=int(dlinex[2]):
				try:
					out_set.add("%s\t%s" % (dline, rlinex[3]))
				except: 
					out_set.add("%s" % dline)
	
#sort_dmc_set = sorted(out_set)
for i in out_set:
	outfile.write("\n%s " % i.strip()) 
dmc_line.close()
roi_line.close()
outfile.close()
