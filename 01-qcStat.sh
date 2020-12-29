#!/bin/bash
##Perform Q20/Q30 percentage calculation > Generate QC Summary

js() {
    jq ".$1" $config | sed 's/"//g'
}

var() {
    if [[ $1 == "subsample" ]]
    then
        echo "Subsample run"
        R1=${sample}-R1.sub
        R2=${sample}-R2.sub
    elif [[ $1 == "run" ]]
    then
        echo "Production run"
        R1=${sample}-R1
        R2=${sample}-R2
    fi
}

qcInfo(){
    if [[ $1 == "before" ]]
    then
        cd Raw_FASTQC
        R1_fastqc=$2_fastqc
        R2_fastqc=$3_fastqc
        unzip '*zip'
    elif [[ $1 == "after" ]]
    then
        cd Trimmed_FASTQC
        R1_fastqc=$2-paired_fastqc
        R2_fastqc=$3-paired_fastqc
        unzip '*zip'
    fi

    #From FastQC report
    total_R1=$(grep "Total Sequences" $R1_fastqc/fastqc_data.txt |cut -f 2)
    total_R2=$(grep "Total Sequences" $R2_fastqc/fastqc_data.txt |cut -f 2)
    length_R1=$(grep "Sequence length" $R1_fastqc/fastqc_data.txt |cut -f 2)
    length_R2=$(grep "Sequence length" $R2_fastqc/fastqc_data.txt |cut -f 2)
    gc_R1=$(grep "%GC" $R1_fastqc/fastqc_data.txt |cut -f 2)
    gc_R2=$(grep "%GC" $R2_fastqc/fastqc_data.txt |cut -f 2)
    qua_R1=$(grep "Per base sequence quality" $R1_fastqc/summary.txt |cut -f 1)
    qua_R2=$(grep "Per base sequence quality" $R2_fastqc/summary.txt | cut -f 1)
    adpt_R1=$(grep "Adapter Content" $R1_fastqc/summary.txt | cut -f 1)
    adpt_R2=$(grep "Adapter Content" $R2_fastqc/summary.txt | cut -f 1)

    #Fromm python q30P script report
    total_R1base=$(grep "total bases" ${R1}_q30P.txt | cut -d " " -f 3 |sed 's/)//g')
    total_R2base=$(grep "total bases" ${R2}_q30P.txt | cut -d " " -f 3 |sed 's/)//g')
    R1_q20base=$(grep "q20 bases" ${R1}_q30P.txt | cut -d " " -f 3 |sed 's/)//g')
    R2_q20base=$(grep "q20 bases" ${R2}_q30P.txt | cut -d " " -f 3 |sed 's/)//g')
    R1_q30base=$(grep "q30 bases" ${R1}_q30P.txt | cut -d " " -f 3 |sed 's/)//g')
    R2_q30base=$(grep "q30 bases" ${R2}_q30P.txt | cut -d " " -f 3 |sed 's/)//g')
    R1_q20P=$(grep "q20 percents" ${R1}_q30P.txt | cut -d " " -f 3 |sed 's/)//g' |xargs printf "%0.2f")
    R2_q20P=$(grep "q20 percents" ${R2}_q30P.txt | cut -d " " -f 3 |sed 's/)//g' |xargs printf "%0.2f")
    R1_q30P=$(grep "q30 percents" ${R1}_q30P.txt | cut -d " " -f 3 |sed 's/)//g' |xargs printf "%0.2f")
    R2_q30P=$(grep "q30 percents" ${R2}_q30P.txt | cut -d " " -f 3 |sed 's/)//g' |xargs printf "%0.2f")

    cd ..
}


config="analysis-profile.json"
samplenames="$(js samplenames)"
mode="$(js mode)"
q30P="$(js q30P)"

for sample in $samplenames
do
    echo "$sample"
    cd $sample
    var $mode
    echo "$R1 $R2"
    python $q30P $R1.fastq.gz >Raw_FASTQC/${R1}_q30P.txt &
    python $q30P $R2.fastq.gz >Raw_FASTQC/${R2}_q30P.txt &
    cd Trimmed
    python $q30P $R1-paired.fastq >../Trimmed_FASTQC/${R1}_q30P.txt &
    python $q30P $R2-paired.fastq >../Trimmed_FASTQC/${R2}_q30P.txt &

    cd ../../

done

while [[ $(ps -u $(whoami) | grep python | wc -l) -ge 1 ]]
do
    sleep 30
done

printf "Sample Name\tTotal R1\tTotal R2\tTotal R1 Bases\tTotal R2 Bases\tR1 Length\tR2 Length\tR1 GC%%\tR2 GC%%\tR1 Quality Assessment\t\
R2 Quality Assessment\tR1 Adapter Assessment\tR2 Adapter Assessment\tR1 Q20 Bases\tR2 Q20 Bases\tR1 Q30 Bases\tR2 Q30 Bases\t\
Total R1\tTotal R2\tTotal R1 Bases\tTotal R2 Bases\tR1 Length\tR2 Length\tR1 GC%%\tR2 GC%%\tR1 Quality Assessment\tR2 Quality Assessment\t\
R1 Adapter Assessment\tR2 Adapter Assessment\tR1 Q20 Bases\tR2 Q20 Bases\tR1 Q30 Bases\tR2 Q30 Bases\n">QC_Stat.txt

echo "Empty QC statistics file created in $(pwd)"

for sample in $samplenames
do
    cd $sample
    var $mode
    ln -s ../QC_Stat.txt .

    qcInfo before $R1 $R2
    printf "\n$sample\t$total_R1\t$total_R2\t$total_R1base\t$total_R2base\t$length_R1\t$length_R2\t$gc_R1\t$gc_R2\t$qua_R1\t$qua_R2\t$adpt_R1\t$adpt_R2\t\
$R1_q20base($R1_q20P%%)\t$R2_q20base($R2_q20P%%)\t$R1_q30base($R1_q30P%%)\t$R2_q30base($R2_q30P%%)">>QC_Stat.txt
    qcInfo after $R1 $R2
    printf "\t$total_R1\t$total_R2\t$total_R1base\t$total_R2base\t$length_R1\t$length_R2\t$gc_R1\t$gc_R2\t$qua_R1\t$qua_R2\t$adpt_R1\t$adpt_R2\t\
$R1_q20base($R1_q20P%%)\t$R2_q20base($R2_q20P%%)\t$R1_q30base($R1_q30P%%)\t$R2_q30base($R2_q30P%%)">>QC_Stat.txt

    rm QC_Stat.txt
    rm -r Raw_FASTQC/*_fastqc
    rm -r Trimmed_FASTQC/*_fastqc

    echo "QC statistics for $sample printed"
    cd ../

done

echo "Complete QC statistics file generated, please find it in $(pwd) "

cut -f 1-3,16-17 QC_Stat.txt | sed 's|\t|,|g' >raw_readsSummary.csv
cut -f 1,18-19,32-33 QC_Stat.txt | sed 's|\t|,|g' >trim_readsSummary.csv

echo "QC summary file (before trimming and after trimming) generated, please find it in $(pwd) "
