#!/bin/bash
##Perform DGE with Cuffdiff

js() {
    jq ".$1" $config | sed 's/"//g'
}


config="analysis-profile.json"
samplenames="$(js samplenames)"
comparison_groups="$(js comparison_groups)"

#for comparison in $comparison_groups;do
for comparison in EDvsT0 EDvsT8;do

    printf "Query name, Query Length, Target name, Target Length, Energy, BasePair\n" > lncRNA-mRNA/sele/$comparison.txt
    sed '1,2d' lncRNA-mRNA/results/$comparison.txt | cut -f1 -d"," --complement | sort -t"," -n -k5 -T /data/siting/tmp | awk -F ',' '!x[$1 FS $3]++' | sort -t"," -k1 -T /data/siting/tmp >> lncRNA-mRNA/sele/$comparison.txt &
    
done
