#!/usr/bin/env bash

module load samtools
module load bcftools


export tumor=$1
export normal=$2
export genome=$3

#  ../bams/WGS/B01P0095AAA03_tumour.bam  ../bams/WGS/B01P0095_ABC03_normal.bam

TARGET="interval-exome.target-267.bed"
ACCESS="hg38_access2k.bed"
OUT=`echo $tumor | sed s/.bam//g`
OUT=`basename $OUT`

echo "$OUT $ACCESS $TARGET"

`date`

python3 cnvkit.py batch $1 --normal $2 \
    --fasta $genome --access data/$ACCESS  \
    -t data/$TARGET \
    --output-reference my_reference.cnn --output-dir results/$OUT \

 
 `date`
