#!/bin/bash
module add bowtie2/2.4.1
module list

sample_list=./sample_data/sample_list.csv
te_file=./sample_data/TEseq_seed.csv
reference=./data/GRCm38_68.fa
bowtie_dir=./bowtie
species=mouse
C=25
T=80
K=25

echo "Step 01"
python ./scripts/Step01_FindContext.py $sample_list $te_file $C $T $K
echo "Step 02"
python ./scripts/Step02_MapContext.py $sample_list $reference $bowtie_dir $species $C
echo "Step 03"
python ./scripts/Step03_MergeSides.py $sample_list $te_file
echo "Step 04"
python ./scripts/Step04_FindZygosity.py $sample_list $te_file $C
