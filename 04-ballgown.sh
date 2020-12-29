#!/bin/bash
##Perform DGE with Cuffdiff

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
samplenames="$(js samplenames)"
parallel_processes="$(js parallel_processes)"
samtools="$(js samtools)"
comparison_groups="$(js comparison_group)"
m_dir="$(js exp_stie_dir_m)"
lnc_dir="$(js exp_stie_dir_lnc)"
pattern="$(js pattern)"
measure="$(js measure)"
m_file="$(js m_file)"
lnc_file="$(js lnc_file)"
main_dir="$(js main_dir)"
dgedir="$(js dgedir)"
RNA_rda="$(js RNA_rda)"
num="$(js num)" # to filter and keep the genes expressed at a reasonable level in at least 20% per treatment condition
org="$(js annotation_lib)"
kegg_org="$(js organism)"
html_tail="$(js tail)"
html_header="$(js header)"
map="$(js map)"

# mkdir ExpressionStieRNA

# mkdir ExpressionStielncRNA

# for sample in $samplenames
# do
#     ln -s $main_dir/$sample/ExpressionStieRNA ExpressionStieRNA/sample$sample

#     ln -s $main_dir/$sample/ExpressionStielncRNA ExpressionStielncRNA/sample$sample
# done

# Rscript bg-load.R $m_dir $pattern $measure $m_file

# Rscript bg-load.R $lnc_dir $pattern $measure $lnc_file
# ids=($(ls ExpressionStieRNA))

# find */ExpressionStielncRNA/ -name "*.gtf" -type f -printf "%p \n" > sample_lst_lncRNA.txt
# for id in ${ids[@]};do
#     sed -i "s/^${id#sample}/${id}\t${id#sample}/" sample_lst_lncRNA.txt
# done

# find */ExpressionStieRNA/ -name "*.gtf" -type f -printf "%p \n" > sample_lst_mRNA.txt

# for id in ${ids[@]};do
#     sed -i "s/^${id#sample}/${id}\t${id#sample}/" sample_lst_mRNA.txt
# done


# for rda in $RNA_rda;do
#     python prepDE.py -i sample_lst_lncRNA.txt -g ${rda%.*}/${rda%.*}_gene_count_matrix.csv -t ${rda%.*}/${rda%.*}_transcript_count_matrix.csv
# done

    
get_group() {
    groups="$(jq ".groups | .$1" $config | sed 's/"//g')"
    n=1
    label=$(echo $groups | sed 's| |,|g')
    for group in $groups
    do
	eval g$n=$group
	n=$(echo "$n+1" | bc )
done
}

# for rda in $RNA_rda
# #for rda in bg_m.rda
# do
# for comparison in "C0vsC8" "T8vsT0" 
# do
#     get_group $comparison
#     echo $g1
#     echo $g2
#     Rscript bg-expr-reciprocal-ts.R $g1 $g2 $dgedir $rda $num
# done
# done

for rda in $RNA_rda
do
for comparison in $comparison_groups
do
    get_group $comparison
#    get_group "indvsnor"
    echo $g1
    echo $g2
    Rscript bg-expr.R $g1 $g2 $rda $num $org $kegg_org $html_tail $html_header $map
done
done

# for type in lnc m;do
# for comparison in $comparison_groups;do
# #    sed -r 's/^(>\S+)\s.*/\1/' MergedGTF${type}RNA_Stie/*.fa > MergedGTF${type}RNA_Stie/merged.fa
#    cat DGE/ID_bg_${type}_${comparison}_transcript_sig.csv | xargs $samtools faidx MergedGTF${type}RNA_Stie/merged.fa  > DGE/FA_${type}RNA_${comparison}_transcript_sig.fa
# done
# done
