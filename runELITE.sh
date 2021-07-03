#!/bin/bash
module add bowtie2/2.4.1
module list

sample_list=./sample_data/sample_list.csv
te_file=./sample_data/TEseq_new_seeds.csv
reference=/pine/scr/a/n/anwica/reference/GRCm38_68.fa
bowtie_dir=./bowtie
species=mouse
C=25
T=80
K=25

python ./scripts/Step01_FindContext.py $sample_list $te_file $C $T $K
python ./scripts/Step02_MapContext.py $sample_list $reference $bowtie_dir $species $C
python ./scripts/Step03_MergeSides.py $sample_list $te_file
python ./scripts/Step04_FindZygosity.py $sample_list $te_file $C
