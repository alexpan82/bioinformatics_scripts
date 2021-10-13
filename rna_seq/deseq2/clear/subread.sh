#!/bin/bash

echo "Input Sample Species:"
read type

# lcRNA uses -s 0
# bulk RNA uses -s 2
/users/PAS1088/osu0230/dev/subread-1.5.1-source/bin/featureCounts -s 0 -T 12 -p -B -C --primary -O\
	-a /fs/project/PAS0472/toolbelt/references/$type/annotation.gtf -o counts.csv */align.sorted.bam

cat counts.csv | fgrep -v '#' | sed 's/\/align.sorted.bam//g' | cut -f 1,7- > counts.tsv
