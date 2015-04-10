#!/bin/bash
## warning: each parition of the data needs to be run in its own analysis folder
set -v
set -e

data=$1
## tag is just to identify which cutoff was used to split the data
tag=allSNPs 

## cutOff
#cutOff=( 1.15 1.20 1.25 1.30 1.35 1.40 1.50 1.60 1.70 1.80 2.00 2.25 2.50 3.00 3.50 10.00 )
#cutOff=( 3.00 10.00 )
#cutOff=( 1.15 35.00 )
#cutOff=( 0.001 35.0 )
cutOff=( 1.5 35.0 )

## log2(cutOff)
#lcutOff=( 0.2016339 0.2630344 0.3219281 0.3785116 0.4329594 0.4854268 0.5849625 0.6780719 0.7655347 0.8479969 1.0000000 1.1699250 1.3219281 1.5849625 1.8073549 3.321928 )
#lcutOff=( 1.5849625 3.321928 )
#lcutOff=( 0.2016339 5.129283 )
# yes, below we will capture all of the data, the range is actually c(-5, 5), but whateva
#lcutOff=( -9.965784 5.129283 )
lcutOff=( 0.5849625 5.129283 )

for ii in {0..0}; do
#for ii in {0..14}; do
    lb=${lcutOff[${ii}]};
    ub=${lcutOff[${ii}+1]};
    echo Processing lower: $lb upper: $ub ii: $ii;
    
    less ${data} | grep -v "snp" | awk -v lb=$lb -v ub=$ub '($8>=lb && $8<ub) || ($8<=-lb && $8>-ub)' | awk '{print $1,$2,$3,$4"\n"$1,$5,$6,$7}' | tr " " "\t" > ../analysis/allPlates_${tag}.txt

    echo Processing complete

done;
