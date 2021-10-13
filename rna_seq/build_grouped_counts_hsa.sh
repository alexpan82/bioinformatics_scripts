bash subread.sh  < <( echo hsa )

ln -s counts.tsv counts.rep.tsv
#cat counts.tsv | python ~osu8725/tools/toolbelt/scripts/replaceGencode.py | sed 's/-[0-9][0-9][0-9]//g' > counts.rep.tsv

time python ~osu8725/tools/toolbelt/scripts/group_for_impute.py ../samples.txt 0 > clear.counts.tsv


