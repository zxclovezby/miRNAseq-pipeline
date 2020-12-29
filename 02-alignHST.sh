#!/bin/bash
##Perform reads alignment with HISAT2

js() {
    jq ".$1" $config | sed 's/"//g'
}

delay() {
    while [[ $(ps -f -u $(whoami) | grep $1 | wc -l) -ge $parallel_processes ]]
    do
        sleep 10
    done
}

hisat() {
    $hisat/hisat2 --mm -p $hisat_threads -x $1 --dta-cufflinks -1 $2 -2 $3 | $samtools view -@ $samtools_threads -Sbh -| $samtools sort -O bam -@ $samtools_threads -T $sample - >$sample.bam &
}


config="analysis-profile.json"
mode="$(js mode)"
samplenames="$(js samplenames)"
parallel_processes="$(js parallel_processes)"
hisat="$(js hisat)"
hisat_threads="$(js hisat_threads)"
hisat_index="$(js hisat_index)"
samtools="$(js samtools)"
samtools_threads="$(js samtools_threads)"

for sample in $samplenames
do
    echo "$sample $mode"
    cd $sample

    if [[ $mode == "subsample" ]]
    then
	echo "Read alignment for $sample subsample at $(date)"
	R1=${sample}-R1.sub
	R2=${sample}-R2.sub
    elif [[ $mode == "run" ]]
    then
	echo "Read alignment for $sample at $(date)"
	R1=${sample}-R1
	R2=${sample}-R2
    fi

    echo "$R1 $R2"

    delay hisat2
    mkdir AlignmentHST
    cd AlignmentHST
    ln -s ../Trimmed/*-paired.fastq .
    hisat $hisat_index $R1-paired.fastq $R2-paired.fastq
    cd ../../

done

while [[ $(ps -u $(whoami) | grep samtools | wc -l) -ge 1 ]]
do
    sleep 10
done

echo "Read alignment completed at $(date)"

for sample in $samplenames
do

    if [[ $mode == "subsample" ]]
    then
        R1=${sample}-R1.sub
        R2=${sample}-R2.sub
    elif [[ $mode == "run" ]]
    then
        R1=${sample}-R1
        R2=${sample}-R2
    fi

    cd $sample/AlignmentHST
    rm $R1-paired.fastq $R2-paired.fastq

    delay samtools
    $samtools flagstat $sample.bam > ${sample}_algnStat.txt &
    cd ../../

done

while [[ $(ps -u $(whoami) | grep samtools | wc -l) -ge 1 ]]
do
    sleep 10
done

echo "Read alignment statistics generated at $(date)"
