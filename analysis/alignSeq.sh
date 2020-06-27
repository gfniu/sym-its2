
#!/bin/bash

SAVEIFS=$IFS
IFS=$(echo -en "\n\b")

# --------------------
# 1 Haplotype alignment
# --------------------
:<< 'COMMENT'
INPUT_FOLDER='/Users/niuguigui/Documents/GitHub/sym-its2/data/5 Haplotype/sequences'
OUTPUT_FOLDER='/Users/niuguigui/Documents/GitHub/sym-its2/data/5 Haplotype/mafft_alignment'
for READ in `ls $INPUT_FOLDER`;
do
     if [[ $READ == *.fasta ]];
     then
           echo "$READ"
           mafft --globalpair --maxiterate 1000 --reorder  ''$INPUT_FOLDER'/'$READ'' > ''$OUTPUT_FOLDER'/'${READ%%.*}'.mafft.aln.fasta'
     fi
done
COMMENT

# --------------------------
# 2 Genetic distance alignment
# --------------------------
:<< 'COMMENT'
INPUT_PATH='/Users/niuguigui/Documents/GitHub/sym-its2/data/7 Genetic distance/Sym_ITS2_core_583.fasta'
OUTPUT_PATH='/Users/niuguigui/Documents/GitHub/sym-its2/data/7 Genetic distance/Sym_ITS2_core_583_aln.fasta'

mafft --globalpair --maxiterate 1000 --reorder $INPUT_PATH> $OUTPUT_PATH
COMMENT

# --------------------
# 3 Phylogeny alignment
# --------------------
:<< 'COMMENT'
INPUT_PATH='/Users/niuguigui/Documents/GitHub/sym-its2/data/6 Phylogeny/Sym_ITS2_extended.fasta'
OUTPUT_PATH='/Users/niuguigui/Documents/GitHub/sym-its2/data/6 Phylogeny/Sym_ITS2_extended_aln.fasta'

mafft --globalpair --maxiterate 1000 --reorder $INPUT_PATH> $OUTPUT_PATH
COMMENT

IFS=$SAVEIFS
