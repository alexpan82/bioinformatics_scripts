# Written on 3/11/20
# Corrected. DO NOT USE V1 (count_ratioCutoff.py)

# 1) run methylkit on all uniq aligned reads and subsampled reads -> *CpG.txt
# 2) run correlation_cov.sh and format.sh on all uniq aligned reads .CpG.txt
# 3) Run this script on output from 2) and emseq and rrbs subsampled.CpG.txt
# Usage: python count_ratioCutoff.py E_R_cov_gauss_norm_100_ignore01.smooth.txt emseq_CpG.txt rrbs_CpG.txt

# Algorithm
# 1) E_R_cov_gauss*.txt contains the ratio vs cov graphs
# 2) For every ratio cutoff (from 0 to 0.2)
	# a) Choose the first point (from 0 cov) that goes PAST cutoff
	# b) Count ALL CpGs (either from emseq_CpG.txt or rrbs_CpG.txt) whose cov >= that point
# 3) Plot # CpGs vs cov cutoff

from collections import defaultdict
from sys import argv


# Loop thru 1st file and store coverage and ratio info in dic
# I will also round ratio to 2 sig figs
# covratio[0][cov] = ratio

covratio = defaultdict(lambda: defaultdict(lambda: 0))

with open(argv[1], 'r') as f:
	for n, line in enumerate(f):
		line = line.strip().split("\t")
		cov0 = int(line[0])
		cov1 = int(line[2])
		ratio0 = round(float(line[1]), 3)
		ratio1 = round(float(line[3]), 3)

		covratio[0][cov0] = ratio0
		covratio[1][cov1] = ratio1

# Hold #CpGs with a certain cov
# covcount[0][cov] = count
covcount = defaultdict(lambda: defaultdict(lambda: 0))

# Loop thru subsampled files
for n, subsampled in enumerate(argv[2:]):
	with open(subsampled, 'r') as f:
		for m, line in enumerate(f):
			line = line.strip().split("\t")
			try:
				cov = int(line[4])
				covcount[n][cov] += 1
			except:
				pass

# The reason we only look at ratios 0-0.2 is because very noisey past 0.2
# I will hard-code this list
ratiocutoff = list(map(lambda x: x/1000, range(0,201)))

# Hold info regarding #CpGs with cov >= than that of the cov associated with a certain ratio 
	# (starting from left)
# ratiocpgsum[0][ratio] = sum of cpgs with cov >= covratio[0][cov]
ratiocpgsum = defaultdict(lambda: defaultdict(lambda: 0))

# Holds the leftmost cutoff at a certain ratio cutoff
# covcutoff[0][ratio] = cov
covcutoff = defaultdict(lambda: defaultdict(lambda: 0))

# Loop thru ratiocutoff
for n in sorted(covcount.keys()):
	# This is a sorted list of coverage counts
	# for example: covlist[0] = # of CpGs w/ cov=1
	covlist = [covcount[n][key] for key in sorted(covcount[n].keys())]

	for r in ratiocutoff:
		for cov,numcpg in enumerate(covlist):
			if covratio[n][cov+1] <= r:
				pass
			else:
				covcutoff[n][r] = cov+1
				ratiocpgsum[n][r] += sum(covlist[cov:])
				break

# print intersection point
for n, r in enumerate(sorted(set(ratiocpgsum[0].keys()) | set(ratiocpgsum[1].keys()))):
	if ratiocpgsum[0][r] <= ratiocpgsum[1][r]:
		print('Approx intersection: %s' % r)
		break

# Print results
print('ratio\temseqcov\temseqcount\trrbscov\trrbscount')
for r in sorted(set(ratiocpgsum[0].keys()) | set(ratiocpgsum[1].keys())):
	print('%s\t%s\t%s\t%s\t%s' % (r, covcutoff[0][r], ratiocpgsum[0][r], covcutoff[1][r], ratiocpgsum[1][r]))













