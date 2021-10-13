# authored by Alex
# parses 10X using Change-O
#
# add --partial to MakeDb.py line to allow partial alignments
# (see Change-O documentation for more info)

source activate ~osu8725/anaconda3 #for access to changeo parsing scripts

seq=$1 # filtered_contig.fasta from cellranger vdj
annotations=$2 # filtered_contig_annotations.csv from cellranger vdj

AssignGenes.py igblast -s ${seq} -b /users/PAS0472/osu8725/tools/igblast \
   --organism human --loci ig --format blast \
   --vdb human_gl_V --ddb human_gl_D --jdb human_gl_J

python ~/tools/MakeDb.py igblast -i ${seq%.fasta}_igblast.fmt7 -s ${seq} \
   -r ~osu8725/tools/igblast/database/IMGT/IG[VDJ].fasta --10x ${annotations} \
   --partial --extended
