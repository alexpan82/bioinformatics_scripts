# Prints all cellbarcodes in 10X bam file
# Ex: python split_10x_cellbarcode.py 10x.bam

import pysam
import gzip
import os
from sys import argv
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--bam","-b", type = str, dest = "bam", help = "Input cellranger count outputted bam file")
args = parser.parse_args()

# Read valid cell barcodes to set
cellbarcodes = set()
barcodeDict = defaultdict(lambda: set())

# TODO: Get name of the read (col 0). line.query_name gets the name
with pysam.AlignmentFile(args.bam, "rb") as f:
	tenxheader = f.header
	for n, line in enumerate(f):
		try:
			chrm = "chr" + str(line.reference_name)
			tags = line.get_tags()
			for t in tags:
				if t[0] == 'CB':
					CBcall = t[1]
					cellbarcodes.add(CBcall)
					break

		except:
			pass

with pysam.AlignmentFile(args.bam, "rb") as f:
	for n, line in enumerate(f):
		try:
			chrm = "chr" + str(line.reference_name)
                        tags = line.get_tags()
                        for t in tags:
                                if t[0] == 'CB':
                                        barcodeDict[t[1]].add(line.query_name)
					break
		except:
			pass


with open ("all_reads.txt", "w") as w:
	for barcode in barcodeDict.keys():
		w.write(barcode + "\t" + str(len(barcodeDict[barcode])))
		w.write("\n")


