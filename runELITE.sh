#!/bin/bash

sample_list=./sample_data/sample_list.csv
te_file=./sample_data/TEseq.csv
reference=./sample_data/Reference.fa
bowtie_dir=/pine/scr/a/n/anwica/bowtie_mouse
species=mouse
C=25
T=80
K=25

python ./scripts/Step01_FindContext.py $sample_list $te_file $C $T $K
python ./scripts/Step02_MapContext.py $sample_list $reference $bowtie_dir $species $C
python ./scripts/Step03_MergeSides.py $sample_list $te_file
python ./scripts/Step04_FindZygosity.py $sample_list $te_file $C
