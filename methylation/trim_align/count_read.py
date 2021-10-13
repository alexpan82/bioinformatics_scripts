from sys import argv
import gzip
script, fastq = argv

fastq_line = gzip.open(fastq, 'rb')
#fastq_line.readline()

i = 1
dic = {}
for line in fastq_line:
	#print line
	if i == 2:
		line = line.strip()
		#print line
		try:
			dic[len(line)] += 1
		except:
			dic[len(line)] = 0
			dic[len(line)] += 1
	if i == 4:
		i = 0
	i += 1

for n in sorted(dic.keys()):
	print "%s\t%s" % (n, dic[n])

