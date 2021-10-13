# RRBS Patterns found by sliding window analysis using test of equal proportions
	### For now, script takes 2 .permeth files (#chr\tposition\tcov\tnumCs\t5mc_ratio) that can be outputs from cpg_fill.py. In addition, user must enter probability of success and define what success is.
### Written by Alex Pan

import argparse
import os
import rpy2
import sys
import rpy2.robjects as robjects
def main(args):

	# Parse arguments into variables
	cpg1 = args.cpg1
	cpg2 = args.cpg2
	prob = float(args.prob)
	success = float(args.success)
	wsize = int(args.wsize)
	comb = args.comb
	out = args.out + ".txt"
	#print comb
	# Write a quick error code.
	if success not in [1,2]:
		sys.exit("Please let -s 1 or -s 2")

	if prob >= 1 or prob <=0:
		sys.exit("0 < probability < 1")

	# Read files and store. Shouldnt be a need for coverage filter, as the steps leading up to this (.permeth -> avg_methyl_reps.py -> cpg_fill.py) have already applied a sort of filter
	# Store info into dictionary, with keys of position and values of 5mc ratio

	if comb == "yes":
		cpg1_line = open(cpg1, 'r')
		cpg1_line.readline()
		cpg1_dic_temp = {}
		cpg1_dic = {}
		# Skip header
		itercpg1 = iter(cpg1_line)
		for line in itercpg1:
			line = line.strip().split('\t')
			chrm = int(line[0][3:])
			ratio = float(line[4])
			position = int(line[1])
			# This conditional filters those -0.01 lines
			if ratio >=0:
				cov = float(line[2])
				numc = float(line[3])
				try:
					cpg1_dic_temp[chrm][position] = (cov,numc)
				except:
					cpg1_dic_temp[chrm] = {}
					cpg1_dic_temp[chrm][position] = (cov,numc)

		cpg1_line.close()

		for c in cpg1_dic_temp.keys():
			cpg1_dic[c] = {}
			for p in sorted(cpg1_dic_temp[c].keys()):
				if (p - 1) in cpg1_dic[c].keys():
					total = (cpg1_dic_temp[c][p - 1][0] + cpg1_dic_temp[c][p][0],cpg1_dic_temp[c][p - 1][1] + cpg1_dic_temp[c][p][1])
					cpg1_dic[c][p - 1] = total[1] / total[0]
					#print "%s\t%s\t%s\t%s" % (cpg1_dic_temp[c][p - 1], cpg1_dic_temp[c][p], total, cpg1_dic[c][p - 1])
				if (p - 1) not in cpg1_dic[c].keys():
					cpg1_dic[c][p] = cpg1_dic_temp[c][p][1] / cpg1_dic_temp[c][p][0]	

		cpg2_line = open(cpg2, 'r')
		cpg2_line.readline()
		cpg2_dic_temp = {}
		cpg2_dic = {}
		# Skip header
		itercpg2 = iter(cpg2_line)
		for line in itercpg2:
			line = line.strip().split('\t')
			chrm = int(line[0][3:])
			ratio = float(line[4])
			position = int(line[1])
			# This conditional filters those -0.01 lines
			if ratio >=0:
				cov = float(line[2])
				numc = float(line[3])
				try:
					cpg2_dic_temp[chrm][position] = (cov,numc)
				except:
					cpg2_dic_temp[chrm] = {}
					cpg2_dic_temp[chrm][position] = (cov,numc)

		cpg2_line.close()

		for c in cpg2_dic_temp.keys():
			cpg2_dic[c] = {}
			for p in sorted(cpg2_dic_temp[c].keys()):
				if (p - 1) in cpg2_dic[c].keys():
					total = (cpg2_dic_temp[c][p - 1][0] + cpg2_dic_temp[c][p][0],cpg2_dic_temp[c][p - 1][1] + cpg2_dic_temp[c][p][1])
					cpg2_dic[c][p - 1] = total[1] / total[0]
				if (p - 1) not in cpg2_dic[c].keys():
					cpg2_dic[c][p] = cpg2_dic_temp[c][p][1] / cpg2_dic_temp[c][p][0]	

	if comb == "no":	
		cpg1_line = open(cpg1, 'r')
		cpg1_line.readline()
		cpg1_dic = {}
		# Skip header
		itercpg1 = iter(cpg1_line)
		for line in itercpg1:
			line = line.strip().split('\t')
			chrm = int(line[0][3:])
			ratio = float(line[4])
			position = int(line[1])
			# This conditional filters those -0.01 lines
			if ratio >=0:
				try:
					cpg1_dic[chrm][position] = ratio
				except:
					cpg1_dic[chrm] = {}
					cpg1_dic[chrm][position] = ratio

		cpg1_line.close()

		cpg2_line = open(cpg2, 'r')
		cpg2_line.readline()
		cpg2_dic = {}
		# Skip header
		itercpg2 = iter(cpg2_line)
		for line in itercpg2:
			line = line.strip().split('\t')
			chrm = int(line[0][3:])
			ratio = float(line[4])
			position = int(line[1])
			# This conditional filters those -0.01 lines
			if ratio >=0:
				try:
					cpg2_dic[chrm][position] = ratio
				except:
					cpg2_dic[chrm] = {}
					cpg2_dic[chrm][position] = ratio

		cpg2_line.close()

	### Find difference between 5mC ratios of CpGs that both files share and determine whether this is a "success".

	# Define binary dictionary. Turn position keys into sets and intersect sets. This needs to be done for every #chr so we loop over the chrm keys in the dictionaries as well. For simplicity's sake, we assume both files represent the same # of chromosomes

	outputfile = open(out, 'w')
	# This calls R-function prop.test
	test = robjects.r['prop.test']
	adjust = robjects.r['p.adjust']
	pdic = {} # Dictionary that holds window position:p-value. This is so we can FDR adjust.
	outdic = {} # Holds each output in element
	pvector = []
	outvector = []

	for chrm in sorted(cpg1_dic.keys()):
		pdic[chrm] = {}
		outdic[chrm] = {}
		cpg1_pos = set(cpg1_dic[chrm].keys())
		cpg2_pos = set(cpg2_dic[chrm].keys())
		shared = sorted(cpg1_pos.intersection(cpg2_pos))
		#print sorted(cpg1_pos)
		binary = {}
		for pos in shared:
			if success == 1:
				diff = cpg1_dic[chrm][pos] - cpg2_dic[chrm][pos]
			if success == 2:
				diff = cpg2_dic[chrm][pos] - cpg1_dic[chrm][pos]

			if diff > 0:
				binary[pos] = 1
			if diff <= 0:
				binary[pos] = 0

	### Make windows using this data (wsize bp windows with at least > 5 CpGs). Then plug each window into stat test
		# This temp window stores the same list as window. When loop starts over and window gets a new range, we want to see if the new window has the same CpG positions as the last, so we compare against temp_window. If they are the same, we will pass because we don't want to do a stat test for the same CpG positions.
		temp_window = []
		for i in range(shared[0], shared[-1] + 1):
			window = [x for x in shared if i <= x <= i + wsize]
			#print "%s, %s" %(i,i+wsize)
			#print window
			b_window = []
			if len(window) >= 5 and window != temp_window:
				for pos in window:
					try:
						b_window.append(binary[pos])
					except:
						pass
			# stat test goes here. "x" # successes is count of sufficiently different CpGs "1" and total "n" length of b_window is number of trials
				x = b_window.count(1)
				n = len(b_window)

			# Window cutoff
				if len(b_window) >= 5:
					pdic[chrm][(i,i+wsize)] = test(x, n, p = prob, alternative = "greater")[2][0]
					outdic[chrm][(i,i+wsize)] = "%s:%s-%s\n%s\n%s\n%s" % ("chr"+str(chrm),i,i+wsize,window,b_window,test(x, n, p = prob, alternative = "greater"))
					#print "%s\t%s\t%s" % (i, pvalue[chrm][(i,i+wsize)], test(x, n, p = prob, alternative = "greater")[2][0])

			elif len(window) < 5 or window == temp_window:
				continue
			temp_window = window

	for c in sorted(pdic.keys()):
		for p in sorted(pdic[c].items()):
			#print p
			pvector.append(float(p[1]))
		for o in sorted(outdic[c].items()):
			#print o[1]
			outvector.append(o[1])
	#sort_out = sorted(outdic[chrm].items())
	#print sort_out
	#print robjects.StrVector(pvector)
	qvector = adjust(robjects.StrVector(pvector),"fdr")
	#print qvector
	count = 0
	for i in range(len(pvector)):
		#print sort_out[i]
		#if i < len(): Need to loop over chrms as well
		if qvector[i] <= 0.05:
			outputfile.write("%s\n%s\n\n" % (outvector[i].rstrip(), "q-value = "+str(qvector[i])))
			count += 1
	outputfile.close()

	#print count




if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Generates ROI file from a list of genes')
	parser.add_argument('-c1', dest='cpg1', type=str, help='.permeth file that contains the same CpGs as cpg2 (5 column, tab deliminated)')
	parser.add_argument('-c2', dest='cpg2', type=str, help='.permeth file that contains the same CpGs as cpg1 (5 column, tab deliminated)')
	parser.add_argument('-p', dest='prob', type=str, help='Probability of success (see test of equal proportions)')
	parser.add_argument('-s', dest='success', type=str, help='Define success. Set argument = 1 if cpg1 > cpg2 is success, or set = 2 if cpg1 < cpg2 is success.')
	parser.add_argument('-o', dest='out', type=str, help='output name (.txt is added to string by default)')
	parser.add_argument('-w', dest='wsize', type=str, help='Sliding window size in bp')
	parser.add_argument('-comb', dest='comb', default="no", type=str, help='Combine CpGs that are 1 bp from each other (Type "yes", default = "no")')
	args=parser.parse_args()
	main(args)
