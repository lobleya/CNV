#!/usr/bin/env bash

module load python
module load R

SAMPLE=$1

python cnvkit.py export seg $SAMPLE.cns \
       --enumerate-chroms \
       -o $SAMPLE.seg

# Run PureCN by providing the *.cnr and *.seg files
Rscript PureCN.R \
        --out $SAMPLE \
        --sampleid $SAMPLE \
        --tumor   $SAMPLE.cnr \
        --segfile $SAMPLE.seg \
        --normal_panel $SAMPLE\_mapping_bias_agilent_v6_hg38.rds \
        --vcf $SAMPLE\_mutect.vcf \
        --statsfile $SAMPLE\_mutect_stats.txt \
        --snpblacklist hg38_mask.bed \
        --genome hg38 \
        --funsegmentation none \
        --force --postoptimize --seed 123


#python cnvkit.py call $SAMPLE.cns -o $SAMPLE.call.cns
#python cnvkit.py call $SAMPLE.cns -y -m threshold -t=-1.1,-0.4,0.3,0.7 -o $SAMPLE.call.cns
#python cnvkit.py call $SAMPLE.cns -y -m clonal --purity 0.65 -o $SAMPLE.cns
#python cnvkit.py call $SAMPLE.cns -y -v $SAMPLE.vcf -m clonal --purity 0.7 -o $SAMPLE.call.cns
