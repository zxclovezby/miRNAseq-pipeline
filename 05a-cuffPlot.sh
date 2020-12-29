#!/bin/bash
##Generate Cuffdiff Plots (Density Plot , Box Plot , Volcano Plot , Scatter Plot , Heatmap)

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

cd DGERNA

for comparison in $comparison_groups
do
	delay Rscript

	echo "Start generating Cuffdiff plots for $comparison at $(date)"
	cd $comparison
	mkdir images
	Rscript ../../dgeIMG.R
	cd ../
done

cd ../DGElncRNA

for comparison in $comparison_groups
do
	delay Rscript

	echo "Start generating Cuffdiff plots for $comparison at $(date)"
	cd $comparison
	mkdir images
	Rscript ../../dgeIMG.R
	cd ../
done



While [[ $(ps -u $(whoami) | grep Rscript | wc -l) -ge 1 ]]
do
    sleep 10
done

echo "Cuffdiff plots generated for all comparison at $(date)"

