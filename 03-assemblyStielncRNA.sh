#!/bin/bash
##Perform transcript assembly witn StringTie (Ref based) -> GTF merging -> Extract transfrags sequences (FASTA)

js() {
    jq ".$1" $config | sed 's/"//g'
}

delay() {
    while [[ $(ps -f -u $(whoami) | grep $1 | wc -l) -ge parallel_processes ]]
    do
        sleep 10
    done
}

config="analysis-profile.json"
mode="$(js mode)"
samplenames="$(js samplenames)"
genome_gtf="$(js genome_gtf)"
lnc_gtf="$(js lncRNA_gtf)"
genome_fasta="$(js genome_fasta)"
stringtie="$(js stringtie)"
stringtie_threads="$(js cufflinks_threads)"
gffcompare="$(js gffcompare)"
cufflinks="$(js cufflinks)"
cufflinks_threads="$(js cufflinks_threads)"
parallel_processes="$(js parallel_processes)"

for sample in $samplenames
do
    delay stringtie
    cd $sample
    mkdir AssemblyStielncRNA
    cd AssemblyStielncRNA

    algn_bam=../AlignmentHST/$sample.bam

    echo "Transcript assembly for $sample at $(date)"

    $stringtie -p $stringtie_threads -G $lnc_gtf \
    -o transcripts.gtf -A $sample.exp -l $sample $algn_bam &

    cd ../
    
    mkdir AssemblyStieRNA
    cd AssemblyStieRNA

    algn_bam=../AlignmentHST/$sample.bam

    echo "Transcript assembly for $sample at $(date)"

    $stringtie -p $stringtie_threads -G $genome_gtf \
    -o transcripts.gtf -A $sample.exp -l $sample $algn_bam &

    cd ../../

done

while [[ $(ps -u $(whoami) | grep stringtie | wc -l) -ge 1 ]]
do
    sleep 10
done

echo "Transcript assembly completed for all samples at $(date)"

GTF Merging

echo "Starting GTF merging at $(date)"

mkdir MergedGTFlncRNA_Stie
find */AssemblyStielncRNA/ -name "transcripts.gtf" -type f -printf "%p \n" >transcripts_file-lncRNA.list
$stringtie --merge -p $stringtie_threads -G $lnc_gtf -o MergedGTFlncRNA_Stie/mergedlncRNA.gtf transcripts_file-lncRNA.list
# #$gffcompare -r $genome_gtf -G -o MergedGTFlncRNA_Stie/merged_compare MergedGTFlncRNA_Stie/mergedlncRNA.gtf

 mkdir MergedGTFRNA_Stie
find */AssemblyStieRNA/ -name "transcripts.gtf" -type f -printf "%p \n" >transcripts_file-mRNA.list
$stringtie --merge -p $stringtie_threads -G $genome_gtf -o MergedGTFRNA_Stie/mergedmRNA.gtf transcripts_file-mRNA.list


echo "GTF merging completed at $(date)"


echo "Extracting transcripts sequences with merged GTF file at $(date)"

$cufflinks/gffread -w MergedGTFlncRNA_Stie/transcripts-lncRNA.fa -g $genome_fasta MergedGTFlncRNA_Stie/mergedlncRNA.gtf &

$cufflinks/gffread -w MergedGTFRNA_Stie/transcripts-mRNA.fa -g $genome_fasta MergedGTFRNA_Stie/mergedmRNA.gtf &

while [[ $(ps -u $(whoami) | grep gffread | wc -l) -ge 1 ]]
do
    sleep 30
done

echo "Transcripts extraction completed at $(date)"

echo "Start estimate transcript abundances at $(date)"

for sample in $samplenames
do
    delay stringtie
    cd $sample
    mkdir ExpressionStielncRNA
    cd ExpressionStielncRNA

    algn_bam=../AlignmentHST/$sample.bam 

    echo "Start estimating for $sample at $(date)"

    $stringtie -e -B -p $stringtie_threads -G ../../MergedGTFlncRNA_Stie/mergedlncRNA.gtf \
    -o $sample.gtf -A $sample.exp $algn_bam &

    cd ../

    mkdir ExpressionStieRNA
    cd ExpressionStieRNA

    algn_bam=../AlignmentHST/$sample.bam 

    echo "Start estimating for $sample at $(date)"

    $stringtie -e -B -p $stringtie_threads -G ../../MergedGTFRNA_Stie/mergedmRNA.gtf \
    -o $sample.gtf -A $sample.exp $algn_bam &

    cd ../../

done

echo "Estimation of transcript abundances completed at $(date)"
