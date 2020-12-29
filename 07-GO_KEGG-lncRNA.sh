#!/bin/bash
##Perform gene ontology (GO) and pathway (KEGG) analysis with topGO & respectively.

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
annotation_lib="$(js annotation_lib)"
organism="$(js organism)"
pq_cutoff="$(js pq_cutoff)"

mkdir GOlncRNA
mkdir GOlncRNA/BP GOlncRNA/MF GOlncRNA/CC
cd GOlncRNA

for comparison in $comparison_groups
do
	ln -s ../DGElncRNA_Stie/$comparison/gene_exp.diff ./$comparison
	echo "Start GO analysis for $comparison at $(date)"
	Rscript ../topGO.R $comparison $annotation_lib $pq_cutoff
	rm $comparison
done

cd ..

mkdir PATHWAYlncRNA
cd PATHWAYlncRNA

for comparison in $comparison_groups
do
        ln -s ../DGElncRNA_Stie/$comparison/gene_exp.diff ./$comparison
        echo "Start KEGG analysis for $comparison at $(date)"
        Rscript ../enrichKEGG.R $comparison $annotation_lib $organism $pq_cutoff
        rm $comparison
done

cd ..


while [[ $(ps -u $(whoami) | grep Rscript | wc -l) -ge 1 ]]
do
    sleep 10
done

echo "GO & KEGG analysis completed for all comparison at $(date)"
