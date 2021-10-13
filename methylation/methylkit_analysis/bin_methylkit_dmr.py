# This script requires user input on the exact location of roi files (#chr\tstart\tstop)
# Have these locations in a text file with each location is on a new line
# Bins the DMC/DMR methylKit output by region and distinguishes by direction

from sys import argv

script, dmc, roi = argv
# Read lines from DMC/DMR ROI location text file inputs
dmc_line = open(dmc, 'r')
roi_location = open(roi, 'r')

#dmc_line.readline()
#roi_location.readline()

roi_loc_list = []
dmc_set = set()
dmc_list = []


for line in dmc_line:
    dmc_set.add(line.strip())

for dmc in dmc_set:
    dmc = dmc.split('\t')
    dmc_list.append(dmc)

for line in roi_location:
    roi_loc_list.append(line.strip())

print "\t\thyper\thypo"
# Loop thru for every ROI location in the text file input.
for loc in roi_loc_list:
    roi_file = open(loc, 'r')
    #outfile = open("hi.txt", 'w')

    hyper = 0
    hypo = 0
    # Read lines from ROI file and store location info

    roi_file.readline()
    roi_list = []
    roi_dic = {}
    roi_set = set()
    out_set = set()

    # roi_set will contain #chr\tstart\tstop\tgenesym
    # Change this line if genesym is in different column than expected
    for line in roi_file:
        line = line.strip().split('\t')
        try:
            string = "%s\t%s\t%s" % (line[0], line[1], line[2])
            roi_dic[string] = [line[0], line[1], line[2]]
	    #string = [line[0], line[1], line[2], line[5]]
        except:
            string = "%s\t%s\t%s" % (line[0], line[1], line[2])
            #string = [line[0], line[1], line[2]]
    	roi_set.add(string)

    for roi in roi_set:
        roi_list.append(roi_dic[roi])

    # Intersect DMC/DMR positions with ROI locations in roi_set. Count these intersections.

    for rline in roi_list:
        for dline in dmc_list:
            if dline[0] == rline[0]:
                if int(rline[1])<=int(dline[1])<=int(rline[2]) or int(dline[1])<=int(rline[1])<=int(dline[2]):
                    out_set.add("%s" % ('\t'.join(dline)))
                    #print "%s" % ('\t'.join(dline))

    for i in out_set:
        i = i.strip().split('\t')
        try:
           if float(i[6]) > 0:
               hyper += 1
           elif float(i[6]) < 0:
                hypo += 1
        except:
            hyper += 1
            hypo += 1
            #pass

    print "%s\t%s\t%s" % (loc, hyper, hypo)
    roi_file.close()
    del roi_set
    del roi_list
    del roi_dic
