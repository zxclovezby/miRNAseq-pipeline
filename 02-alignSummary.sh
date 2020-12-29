#%/bin/bash
##Extract HISAT2 alignment summary statistics

js() {
    jq ".$1" $config | sed 's/"//g'
}

config="analysis-profile.json"
samplenames="$(js samplenames)"

printf "Sample Name,Number of Trimmed Reads (PE),Number of Reads Mapped (PE),Mapping Rate (%%)\n" >algSummary.csv
for sample in $samplenames
do
	cd $sample/AlignmentHST/
	totalMapped=$(cat ${sample}_algnStat.txt | awk 'NR==5' | cut -d " " -f 1)
	secondary=$(cat ${sample}_algnStat.txt | awk 'NR==2' | cut -d " " -f 1)
	totalReads=$(cat ${sample}_algnStat.txt | awk 'NR==6' | cut -d " " -f 1)
        primary=$(echo "$totalMapped-$secondary"|bc)
	mappingRate=$(echo "scale=2;$primary*100/$totalReads"|bc)
	printf "$sample,$totalReads,$primary,$mappingRate\n" >>../../algSummary.csv
	cd ../../ 
done 
