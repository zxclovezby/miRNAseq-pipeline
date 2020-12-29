#!/bin/bash

js() {
    jq ".$1" $config | sed 's/"//g'
}

delay() {
    while [[ $(ps -f -u $(whoami) | grep $1 | wc -l) -ge $parallel_processes ]]
    do
        sleep 10
    done
}

config="analysis-profile.json"
comparison_groups="$(js comparison_group)"
parallel_processes="$(js parallel_processes)"
genome_fasta="$(js genome_fasta)"
RIblast="$(js RIblast)"
blast_thread="$(js blast_thread)"
samtools="$(js samtools)"

# for type in lnc m;do
#     for comparison in $comparison_groups;do
# 	cut -f 9 bg_${type}/Diff_exp/Gene_diff.${comparison}.xls | sed '1d' > bg_${type}/Diff_exp/ID_Gene_diff_${comparison}.txt
# 	# sed -r 's/^(>\S+)\s.*/\1/' MergedGTF${type}RNA_Stie/*.fa > MergedGTF${type}RNA_Stie/merged.fa
# 	cat bg_${type}/Diff_exp/ID_Gene_diff_${comparison}.txt | xargs $samtools faidx MergedGTF${type}RNA_Stie/merged.fa  > bg_${type}/Diff_exp/FA_${comparison}.fa
#     done
# done


if [ ! -d lncRNA-mRNA ]; then
    mkdir -p lncRNA-mRNA
fi


cd lncRNA-mRNA

if [ ! -d lncRNA-mRNA ]; then
    mkdir -p results
fi

for comparison in $comparison_groups;do

#   $RIblast db -i ../bg_m/Diff_exp/FA_${comparison}.fa -o ${comparison}_mRNA.db &
   $RIblast ris -i ../bg_lnc/Diff_exp/FA_${comparison}.fa -o results/${comparison}.txt -d ${comparison}_mRNA.db &
    
done

cd ..

echo "lncRNA target mRNA prediction was done for all comparison at $(date)"
