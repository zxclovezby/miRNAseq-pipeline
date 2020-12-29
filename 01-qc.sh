#!/bin/bash
##Perform FastQC on the raw FASTQ -> Quality & adapter trimming -> FastQC on trimmed FASTQ

js() {
    jq ".$1" $config | sed 's/"//g'
}

delay() {
    while [[ $(ps -f -u $(whoami) | grep $1 | wc -l) -ge $parallel_processes ]]
    do
        sleep 10
    done
}

qc_raw() {

    $fastqc -t $fastqc_threads --outdir Raw_FASTQC $R1.f*q.gz $R2.f*q.gz &

}

trim() {

    java -jar $trimmomatic PE -phred33 -threads $trimmomatic_threads ../$R1.f*q.gz ../$R2.f*q.gz \
    "$R1-paired.fastq" "$R1-unpaired.fastq" "$R2-paired.fastq" "$R2-unpaired.fastq" \
    ILLUMINACLIP:$adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 AVGQUAL:$base_quality_cutoff \
    2> $sample.trimmomatic.Q$base_quality_cutoff.txt

}

qc_trim() {

    $fastqc -t $fastqc_threads --outdir Trimmed_FASTQC Trimmed/$R1-paired.fastq Trimmed/$R2-paired.fastq &

}

config="analysis-profile.json"
fastq_dir="$(js fastq_dir)"
fastqc="$(js fastqc)"
fastqc_threads="$(js fastqc_threads)"
trimmomatic="$(js trimmomatic)"
trimmomatic_threads="$(js trimmomatic_threads)"
base_quality_cutoff="$(js base_quality_cutoff)"
adapters="$(js adapters)"
parallel_processes="$(js parallel_processes)"
samplenames="$(js samplenames)"

mode="$(js mode)"
subsample_size="$(( $(js subsample_size) * 4 ))"


for sample in $samplenames
do
    echo -e "\n $sample $mode"

    if [[ ! -f $sample/$sample-R1.f*q.gz ]]
    then
        [[ ! -d $sample ]] && mkdir $sample
        R1full=$(find $fastq_dir -name "$sample*_1.f*q.gz" -type f -printf "%p " )

        echo "$R1full"
        if [[ $(echo "$R1full" | wc -w) -gt 1 ]]
        then
            echo "Combining $R1full"
            cat $R1full > $sample/$sample-R1.fastq.gz
        else
            echo "Linking $R1full"
            ln -s $R1full $sample/$sample-R1.fastq.gz
        fi
    fi

    if [[ ! -f $sample/$sample-R2.fastq.gz ]]
    then
        [[ ! -d $sample ]] && mkdir $sample
        R2full=$(find $fastq_dir -name "$sample*_2.f*q.gz" -type f -printf "%p ")
        echo "$R2full"
        if [[ $(echo "$R2full" | wc -w) -gt 1 ]]
        then
            echo "Combining $R2full"
            cat $R2full > $sample/$sample-R2.fastq.gz
        else
            echo "Linking $R2full"
            ln -s $R2full $sample/$sample-R2.fastq.gz
        fi
    fi

    cd $sample
    if [[ $mode == "subsample" ]]
    then
	echo "Subsampling $subsample_size reads from $sample"
	zcat ${sample}-R1.fastq.gz | head -n $subsample_size > ${sample}-R1.sub.fastq
	zcat ${sample}-R2.fastq.gz | head -n $subsample_size > ${sample}-R2.sub.fastq
	gzip -f ${sample}-R1.sub.fastq
	gzip -f ${sample}-R2.sub.fastq
	R1=${sample}-R1.sub
	R2=${sample}-R2.sub
	echo "QC for $sample subsample at $(date)"
    elif [[ $mode == "run" ]]
    then
	echo "QC production run for $sample at $(date)"
	R1=${sample}-R1
	R2=${sample}-R2
    fi

    echo "$R1 $R2"

    delay java

    mkdir Raw_FASTQC
    qc_raw

    mkdir Trimmed
    cd Trimmed
    trim
    cd ..

    mkdir Trimmed_FASTQC
    qc_trim

    cd ..

done

while [[ $(ps -u $(whoami) | grep java | wc -l) -ge 1 ]]
do
    sleep 10
done

echo "QC completed at $(date)"
