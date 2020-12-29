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
main_dir="$(js main_dir)"
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


# cd DGE
# for type in bg_m bg_lnc;do
# 	echo "Start generating plots for $type $comparison at $(date)"
# #	mkdir images
# 	cd images
# 	ln -s ../../${type}_all_FPKM.csv ${type}
# 	Rscript ../../Bg-dgeIMG.R $type
# 	rm ${type}
# 	cd ../
# done

cd $main_dir

for type in "bg_m" "bg_lnc";do
    for comparison in $comparison_groups;do
	echo "Start generating plots for $type $comparison at $(date)"
        get_group $comparison
        cd DGE
	Rscript ../Bg-dgeIMG2.R $type $comparison $g1 $g2
        cd $main_dir
    done
done


echo "Plots generated for all comparison at $(date)"
