import csv
import os
import sys
import glob
from collections import defaultdict


###
#   11/17/2014   
#   MODIFIED TO TAKE "permeth.gz" files
#   MODIFIED TO OUTPUT AVERAGE PERCENT METHYLATION   
###


def main():
    date = sys.argv[1]
    path = '.'
    if not date:
        print 'Input date for output file!'
        return
    #cgb_dir = open("cgb_dir_list.txt",'r').read().split('\n')[:-1]
    summary = open(date+'_5x_permeth.feature_summary.txt','w')
    #sample_files = open('full_ovca_samples.list','r').read().split('\n')[:-1]

    cgb_dir = glob.glob( os.path.join('*.roi'))
    dir_list = []
    for dir in cgb_dir:
        dir_list.append(dir.replace('.roi',''))
    
    #summary.write('\t'.join(dir_list)+'\n')
    sample_files = []
    sampcheck = defaultdict(list)
    for i in enumerate(dir_list):
        #print 'adding samples from '+path+'/'+i[1]+'/*.ssroi'
        curlist = glob.glob( os.path.join(path+'/'+i[1]+'/*.ssroi' ))
        sampcheck[i[1]] = curlist
    #check number of samples in each directory
    samplens = []
    for i in sampcheck:
        samplens.append(len(sampcheck[i]))
    sampnum = set(samplens)
    if len(sampnum) > 1:
        print 'Differing number of .ssroi samples!'
        print 'number of samples:',sampnum
        errdir = [ n for n in sampcheck if len(sampcheck[n]) == min(sampnum) ]
        print 'Directory with smallest number of samples: ', map(str,errdir)
        for i in enumerate(dir_list):
            if i[1] in errdir: 
                del dir_list[i[0]]
                print i[1],' Removed from feature summary'

    #fullfiles = glob.glob( os.path.join(path+dir_list[1]+'/','*.ssroi'))
    #for file in fullfiles:
        #sample_files.append(file.split('_qseq')[0])
    """
    ###this section removed for redundancy!###
    headlist = glob.glob( os.path.join(path+'/'+dir_list[0]+'/*.ssroi') )
    headlist = sorted(headlist)
    head = []
    for i in sorted(headlist):
        #s = i.split('_qseq')[0]
        #s = i.split('_')[0]
        #out = s.split('/')[-1]
        s = i.split("/")[-1]
        out = s.split("_")[0]
        head.append(out)
    print "Outputting header using list: ",head
    #print head
    summary.write('\t'.join(head)+'\n')
    """


    for roi in sorted(dir_list):
        print 'Processing ',roi 
        line = roi

        #write header if first ROI
        if roi == sorted(dir_list)[0]:
            sample_list = glob.glob( os.path.join('./'+roi+'/*.ssroi') )
            head = "Region"
            for samp in sorted(sample_list):
                sname = samp.split("/")[-1]
                sname = sname.split("_")[0]
                head = head +"\t"+sname
            print "Outputing header:\n",head,"\n"
            summary.write(head+"\n") 

        sample_list = glob.glob( os.path.join('./'+roi+'/*.ssroi') )
        for samp in sorted(sample_list):
            print "\tSample name: ",str(samp),"\n"
            count = 0
            num = 0
            print samp + ' processing...'
            try:
                try:fin = open(samp, 'r')
                except(IOError):
                    print 'Can not open sample:', samp
                    continue
                reader = csv.reader(fin,delimiter="\t")
                header = reader.next()
                ### Test for header ###
                try:
                    '''
                    #ssnroi/mmroi
                    perf = header.index("ratios")
                    covf = header.index("covs")
                    '''
                    perf = header.index("ratio")
                    covf = header.index("cov")
                except(ValueError):
                    print "No ratio or coverage column found in " + samp
                    print header
                    print "cov,percent column: ",covf,perf
                    line = line + '\tNA'
                    #skip file
                    continue

                ### Process file ###
                #make list for methRatios
                ratlist = []                    
                #make running list of counts of regions that pass cov threshold
                loclist = []                    
                for i in reader:
                    #sample methylation sum
                    sampsum = 0
                    #counter for loci that have CpGs that meet coverage threshold
                    loccnt = 0
#                    i = i.strip("\n").split("\t")
                    if not i: continue
                    '''
                    #ssnroi/mmroi
                    try:
                        cov = i[covf].split(",")
                        per = i[perf].split(",")
                    except(TypeError):
                        print "connot convert one to list: ",i[covf], i[perf]
                    except(IndexError):
                        print i
                    # ssnroi/mmroi only
                    #sum methylation list at each feature
                    for i,p in enumerate(per):
                        #test if CpG passes coverage threshold
                        try:
                            if float(cov[i]) > 5.0:
                                count = count + float(p)
                                num += 1
                        except(ValueError):
                            continue
                    #test for CpGs in list with coverage threshold
                    if num > 0:
                        percent = count / num
                        loccnt += 1
                        sampsum += percent
                    else: 
                        continue
                    '''
                    #ssroi
                    try:
                        cov = float(i[covf])
			print " LOOK HERE: %s" % cov
                        per = float(i[perf])
			print " LOOK HERE 2: %s" % cov
                    except(TypeError):
                        continue

		    except(ValueError):
			continue
			
                    if cov >= 5.0:
			print " Do we reach? "
                        percent = per
                        loccnt += 1
                        sampsum += percent

                #process and output data
		print "loccnt is: %s" % loccnt
                sampmeth = sampsum / loccnt
                ratlist.append(sampmeth)
                loclist.append(loccnt)

                line = line + '\t' + str(sampmeth)


            except(IOError):
                line = line + '\tNA'
                print 'Warning! Sample ',samp,' skipped!'
        summary.write(line + '\n')

    #stats
    print "Samples: ",",".join(sample_list)
    print "IncludedLoci: ",",".join(map(str,loclist))
    print "MethRatios: ",",".join(map(str,ratlist))


main()
