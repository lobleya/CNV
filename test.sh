#!/usr/bin/env bash

#$-cwd
#-l h_rt=20:00:00,h_vmem=80G


module load python

#export tumor="../bams/WGS/B01P0095BAA03_tumour.bam"
#export normal="../bams/WGS/B01P0095_ABC03_normal.bam"
#  ../bams/WGS/B01P0095AAA03_tumour.bam  ../bams/WGS/B01P0095_ABC03_normal.bam
#   ln -s --force $normal Normal.bam
#   ln -s --force $tumor Tumour.bam
   TARGET="interval-exome.target-267.bed"
   ACCESS="hg38_access2k.bed"
   OUT=`echo $tumor | sed s/.bam//g`
   OUT=`basename $OUT`
   echo "$OUT $ACCESS $TARGET"

python3 cnvkit.py batch tumour.bam --normal Normal.bam \
    --fasta /data/BCI-BioInformatics/anna/genomes/fsa/hg38.fa --access data/$ACCESS  \
    -t data/$TARGET \
    --output-reference my_reference.cnn --output-dir results/$OUT \
    --diagram --scatter
 

