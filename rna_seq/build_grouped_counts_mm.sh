bash subread.sh  < <( echo mm )

cat counts.tsv | python ~osu8725/tools/toolbelt/scripts/replaceGencodeMouse.py | sed 's/-[0-9][0-9][0-9]//g' > counts.rep.tsv

#time python ~osu8725/tools/toolbelt/scripts/group_for_impute.py ../samples.txt 0.25 > counts.passed0.75.clear.tsv
time python ~osu8725/tools/toolbelt/scripts/group_for_impute.py ../samples.txt 0 > clear.counts.tsv


