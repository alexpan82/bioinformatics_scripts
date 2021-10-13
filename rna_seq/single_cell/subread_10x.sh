#!/bin/bash

echo "Input Sample Species:"
read type

# lcRNA uses -s 0
# bulk RNA uses -s 2
/users/PAS1088/osu0230/dev/subread-1.5.1-source/bin/featureCounts -s 2 -T 12 -p -B -C --primary -a /fs/project/PAS0472/toolbelt/references/${type}-10x-chr/annotation.gtf -o counts_STAR.csv */align.STAR.sorted.bam

cat counts_STAR.csv | fgrep -v '#' | sed 's/\/align.STAR.sorted.bam//g' | cut -f 1,7- > counts_STAR.tsv
