# authored by Altan Turkoglu (turkoglu.12@osu.edu)
# parses IgBLAST using Change-O
#
# add --partial to MakeDb.py line to allow partial alignments
# (see Change-O documentation for more info)

igblastdir=~osu8725/tools/igblast/

input=$1
output=output/${input%.fasta}

if [ ! -d ${output%%/*} ]; then mkdir ${output%%/*}; fi

source activate ~osu8725/anaconda3 #for access to changeo parsing scripts

~osu8725/tools/igblast/bin/igblastn     -query $input \
                                        -num_threads 40 \
                                        -domain_system imgt \
                                        -ig_seqtype Ig \
                                        -organism human \
                                        -auxiliary_data ${igblastdir}/optional_file/human_gl.aux \
                                        -outfmt '7 std qseq sseq btop' \
                                        -out ${output}.fmt7

MakeDb.py igblast   -i ${output}.fmt7 \
                    -s ${input} \
                    -r ~osu8725/tools/igblast/database/IMGT/IG[VDJ].fasta \
                    --regions --scores --failed

for i in ${output}*.tab; do mv $i ${i%.tab}.tsv; done
