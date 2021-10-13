###### Takes two tab delimenated files containing cpg information (#chr, pos, cov, numCs, ratio). Then identifies the cpgs in each file that the other does not share. Then inputs a negative cov, numCs, ratio value into the missing position. For my purposes, a missing cpg means there was a cov cutoff, so it is not significant. We input a negative value as a marker so that we can still plot both cpg files and compare each cpg

### Need to \sort $i | uniq >> temp && mv temp $i afterwards ... :(

from collections import defaultdict
from sys import argv
script, cpg1, cpg2, cpg3 = argv
#, cpg4, cpg5, cpg6, cpg7, cpg8 = argv

#cpg1_line = open(cpg1, 'r')
#cpg2_line = open(cpg2, 'r')
#cpg3_line = open(cpg3, 'r')
#cpg4_line = open(cpg4, 'r')
#cpg5_line = open(cpg5, 'r')
#cpg6_line = open(cpg6, 'r')
input_list = [cpg1, cpg2, cpg3]
#, cpg4, cpg5, cpg6, cpg7, cpg8]
#outfile = open(out, 'w')

# dictionary with pos:ratio
# dictionary with sample:pos
data = defaultdict( lambda: defaultdict(lambda: "NA\tNA\t-0.01\tNA") )
# set that holds position
cpg_set = set()

for i in input_list:
    with open(i, 'r') as f:
        iterf = iter(f)
        next(iterf)
        for line in iterf:
        #for line in f:
            line = line.strip('\n')
            data[i][line.split('\t')[0] + "\t" + line.split('\t')[1]] = "%s\t%s\t%s\t%s" % (line.split('\t')[2], line.split('\t')[3], line.split('\t')[4], line.split('\t')[5])
            cpg_set.add(line.split('\t')[0] + "\t" + line.split('\t')[1])
        f.close()
print len(cpg_set)
scpg_set = sorted(cpg_set, key=lambda line: float(line.split('\t')[1]))
#sorted(cpg_set)
for n in data:
    print n
    count = 0
    with open("test"+n, 'w') as f:
        f.write("#chr\tpos\tcov\tnumCs\t5mCratio\tgenesym")
        for m in scpg_set:
            #print data[n][m]
            string = m + "\t" + data[n][m]
            #f.write( "\n%s\t%s" % (m, data[n][m]) )
            f.write("\n%s" % string)
            count += 1
        f.close()
    print count

    ### Need to \sort $i | uniq >> temp && mv temp $i afterwards ... :(
