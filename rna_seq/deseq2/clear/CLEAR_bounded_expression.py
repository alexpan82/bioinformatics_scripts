#!/bin/python2

# Authored by Altan Turkoglu (turkoglu.12@osu.edu) and Alex Pan (alexander.pan@osumc.edu)
# forked from group_for_impute.py (Logan Walker)

# This script will replace CLEAR-nonpassing expression as determined 
# by the average of the lowest passing window from CLEAR (length-normalized).

# Currently, pass_fraction is not used (second argument). The length-
# normalized bounded expression is output to stderr. 

# The script currently has two modes (first argument):
#       - stringent
#           *  bounded expressors must be all CLEAR-passing in one sample,
#              but not the other
#       - lenient
#           *  bounded expessors must be all CLEAR-passing in one sample,
#              but not necessarily in the other. 


# Assumes that the following file structure exists:
#       ./*sample*/passed_genes.txt
#       ./*sample*/align.sorted.bed.dat

# Example usage:
#       python /tools/toolbelt/scripts/CLEAR_bounded_expression.py stringent 0 DMSO_2hr VIP152_2hr counts.rep.tsv coldata.txt

#/fs/scratch/PAS0472/osu9900/mRNASeq/steven_s/210430_Pearlly_GSL-PY-2112/output/deseq2/VIP152_2hr_vs_DMSO_2hr

import sys
from collections import defaultdict

#parsing arguments
mode = sys.argv[1]
pass_fraction = float(sys.argv[2]) #only print genes where this fraction of samples is passing or greater. -1 will print all genes
cgroup1 = sys.argv[3]
cgroup2 = sys.argv[4]
colgroup = sys.argv[5]
countsTable = sys.argv[6]
colData = sys.argv[7]

#initializing variables
WINDOW_SIZE = 250

genes = set()
samples = []
sample_map = dict()
bottom_window = defaultdict( lambda: set() )
expression_bound = defaultdict( lambda: float(0) ) #length-normalized average of expression in the bottom WINDOW_SIZE genes
group = defaultdict( lambda: 0 ) #number of times the gene appears in */passed_genes.txt files
group1 = defaultdict( lambda: 0 )
group2 = defaultdict( lambda: 0 )
lengths = defaultdict( lambda: {} ) #gene lengths as determined from .dat file

samples1 = []
samples2 = []

#create samples array based on coldata and comparison groups
t = False
index = 1
for line in open( colData, 'r' ):
    line = line.strip().split('\t')
    if not t:
        index = line.index(colgroup)
        t = True
        continue
    sample = line[0]
    if ( line[index] == cgroup1 ):
        samples1.append(sample)
        samples.append(sample)
    elif ( line[index] == cgroup2 ):
        samples2.append(sample)
        samples.append(sample)

#create set of all genes, a set of the genes for each sample, 
#and count for how many samples the gene passes
for sample in samples:
    sample_genes = set()
    lines = []
    for line in open( '../' + sample + '/passed_genes.txt' ):
        line = line.strip()
        sample_genes.add( line )
        lines.append( line )
    sample_map[ sample ] = sample_genes
    for gene in lines[ -(WINDOW_SIZE): ]:
        bottom_window[ sample ].add( gene )  
    for a in sample_genes:
        genes.add( a )
        group[ a ] += 1 #number of passing
        if sample in samples1:
            group1[ a ] += 1
        elif sample in samples2:
            group2[ a ] += 1

#read in dats file for gene lengths
for sample in samples:
    for line in open( '../' + sample + '/align.sorted.bed.dat', 'r'): #this is a lot of extra i/o and can probably be optimized
        line = line.strip().split("\t")
        length = int(line[2])
        gene = line[3]
        lengths[ sample ][ gene ] = length

header, counts = ( [], {} )

#read in feature counts table
for line in open( countsTable, 'r' ):
    line = line.strip()
    if len( header ) == 0:
        header = line.replace( 'Geneid', '' ).strip().split("\t")
        continue
    else:
        line = line.strip().split('\t')
        gene = line[0]
        line = line[1:]
        if gene in counts.keys(): #check for existing values
            sum1 = 0
            for i in line:
                sum1 += int(i)
            sum2 = 0
            for key in counts[ gene ].keys():
                sum2 += int(counts[ gene ][ key ])
            if sum1 < sum2: continue #skip if the replacement row is lower expressed
        counts[ gene ] = {}
        for i in range( len( header ) ):
            counts[gene][ header[ i ] ] = float( line[ i ] )
            if gene in bottom_window[ header[ i ] ]:
                expression_bound[ header[ i ] ] += float( line[ i ] ) / lengths[ header[ i ] ][ gene ] / WINDOW_SIZE


#print header
out = 'Geneid\t'

for a in header:
    if a in samples:
        out += a
        out += '\t'

print out[ : -1 ]

#print the genes that pass the limit and replace non-CLEAR passing values with expression_bound values
#limit = len( samples ) * pass_fraction #if the number of samples is below this, then it will not be considered passing
limit1 = len(samples1)
limit2 = len(samples2)

if mode == "stringent":
    for gene in sorted( genes ):
        if group1[ gene ] == limit1 and group2[ gene ] == limit2:
            pass
        elif group1[ gene ] == 0 and group2[ gene ] == limit2:
            pass
        elif group1[ gene ] == limit1 and group2[ gene ] == 0:
            pass
        elif group1[ gene ] < limit1 and group2[ gene ] < limit2:
            continue
        elif group1[ gene ] == limit1 and group2[ gene ] > 0:
            continue
        elif group1[ gene ] > 0 and group2[ gene ] == limit2:
            continue
        out = str( gene ) + "\t"
        if gene not in counts:
            continue
        for sample in header:
            if sample not in samples:
                continue
            if (gene in sample_map[ sample ]) and (gene in counts):
                out += str(int(round( counts[ gene ][ sample ] )))
            else: #calculate bounded expression here for nonpassing genes 
                out += str(int(round( expression_bound[ sample ] * lengths[ sample ][ gene ] )))
            out += '\t'
        out = out[ :-1 ]
        #   the below is commented out so that mRNA 'zeroes' that we have artificially made passing in all genes will appear no matter what.
        #	if "\t0\t" not in out: #exclude overlapping errors
        print out
elif mode == "lenient":
    for gene in sorted( genes ):
        if group1[ gene ] == limit1 and group2[ gene ] == limit2:
            pass
        elif group1[ gene ] == 0 and group2[ gene ] == limit2:
            pass
        elif group1[ gene ] == limit1 and group2[ gene ] == 0:
            pass
        elif group1[ gene ] < limit1 and group2[ gene ] < limit2:
            continue
        elif group1[ gene ] == limit1 and group2[ gene ] > 0:
            pass
        elif group1[ gene ] > 0 and group2[ gene ] == limit2:
            pass
        out = str( gene ) + "\t"
        if gene not in counts:
            continue
        for sample in header:
            if sample not in samples:
                continue
            if (gene in sample_map[ sample ]) and (gene in counts):
                out += str(int(round( counts[ gene ][ sample ] )))
            else: #calculate bounded expression here for nonpassing genes 
                out += str(int(round( expression_bound[ sample ] * lengths[ sample ][ gene ] )))
            out += '\t'
        out = out[ :-1 ]
        #   the below is commented out so that mRNA 'zeroes' that we have artificially made passing in all genes will appear no matter what.
        #	if "\t0\t" not in out: #exclude overlapping errors
        print out

#print bounded expressions to stderr
out = ''
for a in header:
    if a in samples:
        out += a
        out += '\t'
sys.stderr.write( out[ : -1 ] )

out = '\n'
for sample in header:
    if sample not in samples:
        continue
    out += str( expression_bound[ sample ] )
    out += '\t'
sys.stderr.write( out[ : -1 ] + '\n' )
