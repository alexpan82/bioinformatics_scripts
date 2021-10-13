from sys import argv
# Takes 2 methylKit dmc files (#chr\tstart\tend\tmeth.dff in the 1st, 2nd, 3rd, and 7th cols respectively) 
# Attempts to determine those dmcs that are affected by the condition and then modified by treatment. In this case, inputting WTvsCF dmcs and CFvsCFEGCG dmcs will tell me those dmcs that are shared between both lists.


script, dmc1, dmc2 = argv#, dmc3 = argv

dmc1_line = open(dmc1, 'r')
dmc2_line = open(dmc2, 'r')
#dmc3_line = open(dmc3, 'r')
out1file = open("WTvCF_up_CFvCFEGCG_up.txt" , 'w')
out2file = open("WTvCF_up_CFvCFEGCG_down.txt",'w')
out3file = open("WTvCF_down_CFvCFEGCG_up.txt",'w')
out4file = open("WTvCF_down_CFvCFEGCG_down.txt",'w')
dmc1_line.readline()
dmc2_line.readline()
#dmc3_line.readline()

### [pos:[cov, numCs, genesym]]

dmc1_dic = {}
dmc2_dic = {}
#dmc3_dic = {}
out1_set = set()
out2_set = set()
out1file.write("#chr\tstart\tend\tstrand\tpvalue\tqvalue\tmeth.diff\tgenesym\n")
out2file.write("#chr\tstart\tend\tstrand\tpvalue\tqvalue\tmeth.diff\tgenesym\n")

for line in dmc1_line:
	cline = line.strip().split('\t')
	dmc1_dic[cline[0] + "\t" + cline[1] + "\t" + cline[2]] = line
for line in dmc2_line:
	cline = line.strip().split('\t')
	dmc2_dic[cline[0] + "\t" + cline[1] + "\t" + cline[2]] = line
#for line in dmc3_line:
#        cline = line.strip().split('\t')
#        dmc3_dic[cline[0] + "\t" + cline[1] + "\t" + cline[2]] = line

keys_1 = set(dmc1_dic.keys())
keys_2 = set(dmc2_dic.keys())
#keys_3 = set(dmc3_dic.keys())

intersection = keys_1 & keys_2# & keys_3

print len(intersection)
for line in intersection:
	#out1_set.add(dmc1_dic[line])
	#out2_set.add(dmc2_dic[line])
	diff1 = float(dmc1_dic[line].strip().split('\t')[6])
	diff2 = float(dmc2_dic[line].strip().split('\t')[6])
	if diff1 > 0 and diff2 > 0:
		out1file.write("%s" % dmc2_dic[line])
	if diff1 > 0 and diff2 < 0:
		out2file.write("%s" % dmc2_dic[line])
	if diff1 < 0 and diff2 > 0:
		out3file.write("%s" % dmc2_dic[line])
	if diff1 < 0 and diff2 < 0:
		out4file.write("%s" % dmc2_dic[line])
	
#for i in out1_set:
	#out1file.write("\n%s " % i.strip()) 
#for i in out2_set:
	#out2file.write("\n%s " % i.strip()) 

