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
comparison_groups="$(js comparison_groups)"
parallel_processes="$(js parallel_processes)"
genome_fasta="$(js genome_fasta)"
RIblast="$(js RIblast)"
blast_thread="$(js blast_thread)"


if [ ! -d lncRNA-mRNA ]; then
    mkdir -p lncRNA-mRNA
fi

cd lncRNA-mRNA
for comparison in $comparison_groups;do
    
#    ls /data/siting/lncRNAseq-longsee-20180208/DGElncRNA_Stie/$comparison/x* | parallel -j $blast_thread "perl $LncTar/LncTar.pl -p 1 -l {} -m /data/siting/lncRNAseq-longsee-20180208/DGE/siD2vsNC/gene_exp_sig.fa -d -0.08 -s T -o /data/siting/lncRNAseq-longsee-20180208/test3/{/}.tsv"
    #    perl $LncTar/LncTar.pl -p 1 -l ../DGElncRNA_Stie/$comparison/gene_exp_sig.fa -m ../DGE/$comparison/gene_exp_sig.fa -d $lowest -s T -o $comparison-coexp.csv
    # $RIblast db -i /data/siting/lncRNAseq-rat-20180405/oldTest/mRNA.fa -o /data/siting/lncRNAseq-rat-20180405/oldTest/mRNA.db
    # $RIblast ris -i /data/siting/lncRNAseq-rat-20180405/oldTest/query.fa -o result -d /data/siting/lncRNAseq-rat-20180405/oldTest/mRNA.db    
#    $RIblast db -i ../DGE/FA_mRNA_${comparison}_transcript_sig.fa -o ${comparison}_mRNA.db &
    $RIblast ris -i ../DGE/FA_lncRNA_${comparison}_transcript_sig.fa -o results/${comparison}.txt -d ${comparison}_mRNA.db &
    
done

# files=($(ls x*.csv))

# printf "Query\tLength Query\tTarget\tLength Target\tdG\tndG\tStart Position Query\tEnd Position Query\tStart Position Target\tEnd Position Target\n" > lncRNA-mRNA.tsv

# for file in "${files[@]}";do
#     n=$(wc -l < $file)
#     for ((i=2; i<=$n; i+=5 ));do
#         sed -n $i,${i}p $file >> lncRNA-mRNA.tsv
#     done
# done

# printf "Query\tLength Query\tTarget\tLength Target\tdG\tndG\tStart Position Query\tEnd Position Query\tStart Position Target\tEnd Position Target\n" > lncRNA-mRNA-Align.tsv

# for file in "${files[@]}";do
#     n=$(wc -l < $file)
#     sed -n 2,${n}p $file >> lncRNA-mRNA-Align.tsv
# done

# rm x*.tsv
# rm temp*.txt

cd ..

echo "lncRNA target mRNA prediction was done for all comparison at $(date)"

