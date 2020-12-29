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

get_bam() {
    groups="$(jq ".groups | .$1" $config | sed 's/"//g')"
    g1_bam=""
    g2_bam=""
    n=1
    label=$(echo $groups | sed 's| |,|g')

    for group in $groups
    do
	eval g$n=$group
	samples="$(jq ".samples | .$group" $config | sed 's/"//g')"

	for sample in $samples
	do
		if [[ $n == 1 ]]
                then
                   g1_bam="$g1_bam ./$sample/AlignmentHST/$sample.bam"
                else
                   g2_bam="$g2_bam ./$sample/AlignmentHST/$sample.bam"
                fi
	done

	n=$(echo "$n+1" | bc )
    done

    g1_bam=$(echo $g1_bam |sed 's| |,|g')
    g2_bam=$(echo $g2_bam |sed 's| |,|g')
    bam_files="$g1_bam $g2_bam"
}

config="analysis-profile.json"
samplenames="$(js samplenames)"
parallel_processes="$(js parallel_processes)"
cufflinks="$(js cufflinks)"
cufflinks_threads="$(js cufflinks_threads)"
comparison_groups="$(js comparison_groups)"

export PATH=$PATH:$cufflinks

mkdir DGElncRNA

mkdir DGERNA

for comparison in $comparison_groups
do
	delay cuffdiff
	get_bam $comparison
	echo "Starting DGE analysis for $comparison $(date)"
	echo -e "$g1\t: $g1_bam"
	echo -e "$g2\t: $g2_bam"
	echo -e "\n"

	$cufflinks/cuffdiff -o ./DGElncRNA/$comparison -u -L $label -p $cufflinks_threads ./MergedGTFlncRNA/merged.gtf $bam_files &

	$cufflinks/cuffdiff -o ./DGERNA/$comparison -u -L $label -p $cufflinks_threads ./MergedGTFRNA/merged.gtf $bam_files &

done

while [[ $(ps -u $(whoami) | grep cuffdiff | wc -l) -ge 1 ]]
do
    sleep 10
done

echo "DGE analysis completed for all comparison at $(date)"

