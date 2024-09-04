#!/bin/bash
for i in *_10kb.cool
do
echo $i
threads=14
file=${i%_10kb.cool}  

#cooler balance $i
cool5=$file$'_50kb.cool'
cooler coarsen -k 5 -o $cool5 $i 

cooler balance $cool5
cooltools eigs-cis $cool5 --phasing-track /data/zuowu/mm10_gc_cov_50kb.tsv -o $file$'_50kb' --bigwig
cooltools insulation $cool5 --window-pixels 20 > $file$'_50kb_insulation_1Mb.tsv'
cooltools expected-cis $cool5 -o $file$'_50kb_cis_expected.tsv' -p $threads 
cooltools expected-trans $cool5 -o $file$'_50kb_trans_expected.tsv' -p $threads 
cooltools saddle $cool5 $file$'_50kb.cis.vecs.tsv' $file$'_50kb_cis_expected.tsv' -o $file$'_50kb_saddle_E1' --fig pdf --qrange 0.02 0.98 

coolfiles=$file$'_1Mb.cool'
cooler coarsen -k 100 -o $coolfiles $i 
cooler balance $coolfiles 
done
