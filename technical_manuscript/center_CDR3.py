#Takes UMI overlapped fastq file and CDR3 sequence
#Ex: python center_CDR3.py raw-C_atleast-2.fastq.gz CCCC

from sys import argv
from collections import defaultdict
import gzip
import io
import string

sequences = []
mod_seq = []
nuc_dic = defaultdict(lambda: defaultdict(lambda: 0))
counter = 0

def revcomp (seq):
	trans = string.maketrans('ATGC','TACG')
	return seq.translate(trans)[::-1]
'''
Returns reverse complement of UMI strand
'''


with open (argv[1], "r") as f:
	count = 0
	longest = 0
	longest_seq = "A"
	for i in f:
		if argv[2] in i:
			sequences.append(i)
			if i.index(argv[2]) > longest:
				longest = i.index(argv[2])
			if len(i) > len(longest_seq):
				longest_seq = i
		elif argv[2] in revcomp(i):
			sequences.append(revcomp(i))
			if revcomp(i).index(argv[2]) > longest:
				longest = revcomp(i).index(argv[2])
			if len(i) > len(longest_seq):
				longest_seq = i

for i in sorted(sequences):
	adjust = longest - i.index(argv[2])
	modified = " "*adjust + i.strip()
	mod_seq.append(modified)


'''
Creates a list containing the sequences with the centered CDR3's
'''
counter = 0

for i in range(0, len(longest_seq)):
	counter += 1
	for s in mod_seq:
		if len(s) > i:
			nuc_dic[counter][s[i]] += 1
		else:
			break
		

'''
counts the number of base pairs for each column
'''


w = open("aligned_UMI.txt", "w+")

for i in ["A","C","T", "G"]:
	w.write(i + ": ")
	for j in sorted(nuc_dic.keys()):
		w.write(str(nuc_dic[j][i]) + ",")
	w.write("\n")	
	
for i in sorted(mod_seq):
	w.write(" "*3 + i + "\n")

w.close()

