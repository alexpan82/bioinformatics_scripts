import sys
from collections import defaultdict

#parsing arguments
samples_list = sys.argv[1]
pass_fraction = float(sys.argv[2]) #only print genes where this fraction of samples is passing or greater. -1 will print all genes

#initializing variables
MISSING_CHAR = 'NA'

genes = set()
samples = []
sample_map = dict()
group = defaultdict( lambda: 0 ) #number of times the gene appears in */passed_genes.txt files

#create samples array
for line in open( samples_list, 'r' ):
	line = line.strip()
	samples.append( line )

#create set of all genes, a set of the genes for each sample, and count for how many samples the gene passes
for sample in samples:
	sample_genes = set()
	for line in open( sample + '/passed_genes.txt' ):
		line = line.strip()
		sample_genes.add( line )
	sample_map[ sample ] = sample_genes
	for a in sample_genes:
		genes.add( a )
		group[ a ] += 1

header, counts = ( [], {} )

#read in feature counts table
for line in open( 'counts.rep.tsv', 'r' ):
    line = line.strip()
    if len( header ) == 0:
        header = line.replace( 'Geneid', '' ).strip().split("\t")
        continue
    else:
        line = line.strip().split('\t')
        gene = line[0]
        line = line[1:]
        counts[ gene ] = {}
        for i in range( len( header ) ):
            counts[gene][ header[ i ] ] = float( line[ i ] )

#print header
out = 'Geneid\t'

for a in header:
    if a in samples:
        out += a
        out += '\t'

print out[ : -1 ]

#print the genes that pass the limit and replace non-CLEAR passing values with MISSING_CHAR
limit = len( samples ) * pass_fraction #if the number of samples is below this, then it will not be considered passing.

for gene in sorted( genes ):
    if group[ gene ] < limit:
        continue
    out = str( gene ) + "\t"
    if gene not in counts:
        continue
    for sample in header:
        if sample not in samples:
            continue
        if (gene in sample_map[ sample ]) and (gene in counts):
            out += str( counts[ gene ][ sample ] )
        else:
            out += MISSING_CHAR
        out += '\t'
    out = out[ :-1 ]
    #   the below is commented out so that mRNA 'zeroes' that we have artificially made passing in all genes will appear no matter what.
    #	if "\t0\t" not in out: #exclude overlapping errors
    print out
