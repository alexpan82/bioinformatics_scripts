#!/bin/python

# Authored by Altan Turkoglu (turkoglu.12@osu.edu)

# This script combines Bismark alignment reports for samples with
# multiple sequencing runs and alignment files.

# This was made with Bismark v0.17.0 alignments reports in mind.
# The input will be a directory containing alignment reports, and
# an output directory must be specified as well.

# Alignment reports must be in the same directory, and of the form:
# SAMPLE_*_bismark_bt2_PE_report.txt,

# NB:
# This will overwrite any files in the output directory !!
# This script will not work without properly named alignment reports !!

import argparse
import os, sys
import fnmatch

# parse arguments
parser = argparse.ArgumentParser(prog='bismark_report_combiner.py', description='This script combines Bismark alignment reports for samples multiple sequencing runs and alignment files. It assumes that input files are named as SAMPLE_*_bismark_bt2_PE_report.txt.')
parser.add_argument('input', metavar='in-dir', help='an input directory containing Bismark alignment reports')
parser.add_argument('output', metavar='out-dir', help='an output directory for Bismark alignment reports to go')
try:
    args = parser.parse_args()
except:
    sys.exit(0)

#create output directory if it does not already exist
if not os.path.exists(args.output):
    os.makedirs(args.output)

# aggregate file and sample list
files = []
samples = []

for filename in os.listdir(args.input):
    if filename.endswith("_bismark_bt2_PE_report.txt"):
        files.append(filename)
        
for i in files:
    samples.append(i.split('_')[0])

samples = sorted(set(samples))

# iterate over all files for each sample
for sample in samples:
    reports = fnmatch.filter(files, sample + '*')
    outfile = args.output + '/' + sample + '_combined' + str(len(reports)) + '_report.txt'
    
    if len(reports) == 1:
        outfile = args.output + '/' + reports[0]
    
    sequence_pairs = 0
    unique_best_hits = 0
    unaligned = 0
    multimapped = 0
    discarded = 0
    
    ctgact = 0
    gactct = 0
    gactga = 0
    ctgaga = 0
    
    total_c = 0
    
    total_m_cpg = 0
    total_m_chg = 0
    total_m_chh = 0
    total_m_u = 0
    
    total_u_cpg = 0
    total_u_chg = 0
    total_u_chh = 0
    total_u_u = 0
    
    time = 0
    
    for report in reports:
        with open((args.input + '/' + report), 'r') as f:
            for line in f:
                if "Sequence pairs analysed in total:" in line:
                    sequence_pairs += int(line.split("\t")[1])
                if "Number of paired-end alignments with a unique best hit:" in line:
                    unique_best_hits += int(line.split("\t")[1])
                if "Sequence pairs with no alignments under any condition:" in line:
                    unaligned += int(line.split("\t")[1])
                if "Sequence pairs did not map uniquely:" in line:
                    multimapped += int(line.split("\t")[1])
                if "Sequence pairs which were discarded because genomic sequence could not be extracted:" in line:
                    discarded += int(line.split("\t")[1])

                if "CT/GA/CT:" in line:
                    ctgact += int(line.split("\t")[1])
                if "GA/CT/CT:" in line:
                    gactct += int(line.split("\t")[1])
                if "GA/CT/GA:" in line:
                    gactga += int(line.split("\t")[1])
                if "CT/GA/GA:" in line:
                    ctgaga += int(line.split("\t")[1])

                if "Total number of C's analysed:" in line:
                    total_c += int(line.split("\t")[1])

                if "Total methylated C's in CpG context:" in line:
                    total_m_cpg += int(line.split("\t")[1])
                if "Total methylated C's in CHG context:" in line:
                    total_m_chg += int(line.split("\t")[1])
                if "Total methylated C's in CHH context:" in line:
                    total_m_chh += int(line.split("\t")[1])
                if "Total methylated C's in Unknown context:" in line:
                    total_m_u += int(line.split("\t")[1])

                if "Total unmethylated C's in CpG context:" in line:
                    total_u_cpg += int(line.split("\t")[1])
                if "Total unmethylated C's in CHG context:" in line:
                    total_u_chg += int(line.split("\t")[1])
                if "Total unmethylated C's in CHH context:" in line:
                    total_u_chh += int(line.split("\t")[1])
                if "Total unmethylated C's in Unknown context:" in line:
                    total_u_u += int(line.split("\t")[1])

                if "Bismark completed in" in line:
                    time += int(line.split(" ")[6].replace("s", "").strip())
                    time += (int(line.split(" ")[5].replace("m", "").strip()) * 60)
                    time += (int(line.split(" ")[4].replace("h", "").strip()) * 60 * 60)
                    time += (int(line.split(" ")[3].replace("d", "").strip()) * 60 * 60 * 24)

    #begin writing file output
    f = open(outfile, 'w')

    f.write("Bismark report for: " + sample + " (version v0.17.0)\n")
    f.write("Bismark was run with Bowtie 2 against the bisulfite genome of /fs/scratch/osu7905/bs-index/GRCh38/noUnplacedContigs/ with the specified options: -q -N 1 --score-min L,-0.6,-0.6 -p 2 --reorder --ignore-quals --no-mixed --no-discordant --dovetail --maxins 1000\n")
    f.write("Option '--non_directional' specified: alignments to all strands were being performed (OT, OB, CTOT, CTOB)\n")
    f.write("\n")
    f.write("Final Alignment report\n")
    f.write("======================\n")
    f.write("Sequence pairs analysed in total:\t" + str(sequence_pairs) + "\n")
    f.write("Number of paired-end alignments with a unique best hit:\t" + str(unique_best_hits) + "\n")
    f.write("Mapping efficiency:\t" + "%.1f%%" % (100. * unique_best_hits/sequence_pairs) + "\n")
    f.write("Sequence pairs with no alignments under any condition:\t" + str(unaligned) + "\n")
    f.write("Sequence pairs did not map uniquely:\t" + str(multimapped) + "\n")
    f.write("Sequence pairs which were discarded because genomic sequence could not be extracted:\t" + str(discarded) + "\n")
    f.write("\n")
    f.write("Number of sequence pairs with unique best (first) alignment came from the bowtie output:\n")
    f.write("CT/GA/CT:\t" + str(ctgact) + "\t((converted) top strand)\n")
    f.write("GA/CT/CT:\t" + str(gactct) + "\t(complementary to (converted) top strand)\n")
    f.write("GA/CT/GA:\t" + str(gactga) + "\t(complementary to (converted) bottom strand)\n")
    f.write("CT/GA/GA:\t" + str(ctgaga) + "\t((converted) bottom strand)\n")
    f.write("\n")
    f.write("Final Cytosine Methylation Report\n")
    f.write("=================================\n")
    f.write("Total number of C's analysed:\t" + str(total_c) + "\n")
    f.write("\n")
    f.write("Total methylated C's in CpG context:\t" + str(total_m_cpg) + "\n")
    f.write("Total methylated C's in CHG context:\t" + str(total_m_chg) + "\n")
    f.write("Total methylated C's in CHH context:\t" + str(total_m_chh) + "\n")
    f.write("Total methylated C's in Unknown context:\t" + str(total_m_u) + "\n")
    f.write("\n")
    f.write("\n")
    f.write("Total unmethylated C's in CpG context:\t" + str(total_u_cpg) + "\n")
    f.write("Total unmethylated C's in CHG context:\t" + str(total_u_chg) + "\n")
    f.write("Total unmethylated C's in CHH context:\t" + str(total_u_chh) + "\n")
    f.write("Total unmethylated C's in Unknown context:\t" + str(total_u_u) + "\n")
    f.write("\n")
    f.write("\n")
    f.write("C methylated in CpG context:\t" + "%.1f%%" % (100. * total_m_cpg/(total_m_cpg + total_u_cpg)) + "\n")
    f.write("C methylated in CHG context:\t" + "%.1f%%" % (100. * total_m_chg/(total_m_chg + total_u_chg)) + "\n")
    f.write("C methylated in CHH context:\t" + "%.1f%%" % (100. * total_m_chh/(total_m_chh + total_u_chh)) + "\n")
    f.write("C methylated in unknown context (CN or CHN):\t" + "%.1f%%" % (100. * total_m_u/(total_m_u + total_u_u)) + "\n")
    f.write("\n")
    f.write("\n")
    f.write("Bismark completed in " + str(time / 86400) + "d " + str((time % 86400) / 3600) + "h " + str((time % 3600) / 60) + "m " + str(time % 60) + "s\n")

    f.close()

print 'Bismark_report_combiner finished running successfully.'
sys.exit(1)