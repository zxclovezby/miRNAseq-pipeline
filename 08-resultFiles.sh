#!/bin/bash
##Link or copy result files into result folder

js() {
    jq ".$1" $config | sed 's/"//g'
}

fastq_plot() {

    #Rscript $main_dir/ngs-qc-plots.R . $1
    unzip -c ${1}_fastqc.zip ${1}_fastqc/fastqc_data.txt > fastqc_data.txt
    python $main_dir/fastqc-split.py fastqc_data.txt

    if [ $(cat length-dist.txt |wc -l) -eq 1 ]
    then
        length1=$(echo "$(cut -f 1 length-dist.txt )-1" |bc -l)
        length2=$(echo "$(cut -f 1 length-dist.txt )+1" |bc -l)
        sed -i "1s|^|$length1\t0\n|g" length-dist.txt
        echo -e "$length2\t0\n" >>length-dist.txt
    fi

    Rscript $main_dir/ngs-qc-plots.R . $1
}

config="analysis-profile.json"
replicate_status="$(js replicate_status)"
samplenames="$(js samplenames)"
main_dir="$(js main_dir)"
comparison_groups="$(js comparison_groups)"

# ##QC files
# for sample in $samplenames
# do
# 	cd $sample/Trimmed
# 	gzip *.fastq &
# 	cd $main_dir
# done

# while [[ $(ps -u $(whoami) | grep gzip | wc -l) -ge 1 ]]
# do
#     sleep 10
# done

# mkdir $main_dir/Report
# mkdir $main_dir/Report/Result
 cd $main_dir/Report/Result
# mkdir QC Alignment Assembly Merge StringTieExp DiffExp GO Pathway lncRNA-mRNA


# cd QC
# mkdir Raw_QC Trimmed Trimmed_QC
# mkdir Raw_QC/Images Trimmed_QC/Images
# for sample in $samplenames
# do
# 	cp $main_dir/$sample/Raw_FASTQC/${sample}-R1_fastqc.* ./Raw_QC
#         cp $main_dir/$sample/Raw_FASTQC/${sample}-R2_fastqc.* ./Raw_QC
#         cp $main_dir/$sample/Trimmed_FASTQC/${sample}-R1-paired_fastqc.* ./Trimmed_QC
#         cp $main_dir/$sample/Trimmed_FASTQC/${sample}-R2-paired_fastqc.* ./Trimmed_QC

#         cd ./Raw_QC
#         fastq_plot $sample-R1
#         fastq_plot $sample-R2
# 	mv *-base-quality.png Images
#         mv *-read-quality.png Images
#         mv *-length-dist.png Images
#         rm base-quality.txt length-dist.txt fastqc_data.txt read-quality.txt Rplots.pdf 
#         cd ..

#         cd ./Trimmed_QC
#         fastq_plot $sample-R1-paired
#         fastq_plot $sample-R2-paired
#         mv *-base-quality.png Images
#         mv *-read-quality.png Images
#         mv *-length-dist.png Images
# 	rm base-quality.txt length-dist.txt fastqc_data.txt read-quality.txt Rplots.pdf
# 	cd ..

#         mkdir ./Trimmed/$sample
#         cd ./Trimmed/$sample
#         ln -s $main_dir/$sample/Trimmed/* .
#         rm *trimmomatic.Q30.txt

#         cd $main_dir/Report/Result/QC
# done
# cd ..

# #Alignment files
# cd Alignment

# for sample in $samplenames
# do
#    ln -s $main_dir/$sample/AlignmentHST/${sample}.bam .
# done
# cd ..

# #Assembly files
# cd Assembly

# for sample in $samplenames
# do
#     ln -s $main_dir/$sample/AssemblyStieRNA/transcripts.gtf ./$sample-mRNA.gtf
#     ln -s $main_dir/$sample/AssemblyStielncRNA/transcripts.gtf ./$sample-lncRNA.gtf
# done
# cd ..

# #Assemblies merged files
# cd Merge

#    ln -s $main_dir/MergedGTFmRNA_Stie/mergedmRNA.gtf ./merged-mRNA.gtf
#    ln -s $main_dir/MergedGTFmRNA_Stie/transcripts-mRNA.fa ./transcripts-mRNA.fa
#    ln -s $main_dir/MergedGTFlncRNA_Stie/mergedlncRNA.gtf ./merged-lncRNA.gtf
#    ln -s $main_dir/MergedGTFlncRNA_Stie/transcripts-lncRNA.fa ./transcripts-lncRNA.fa   
# cd ..

# ## StringTie expression
 cd StringTieExp

# ln -s $main_dir/ExpressionStielncRNA lncRNA
# ln -s $main_dir/ExpressionStieRNA mRNA

# mkdir Images

for types in bg_lnc bg_m
do
    ln -s $main_dir/${types}/Sample_QC/Box*.pdf Images/${types}_Boxplot.pdf
    ln -s $main_dir/${types}/Sample_QC/Box*.png Images/${types}_Boxplot.png
    ln -s $main_dir/${types}/Sample_QC/Cluster*.pdf Images/${types}_Cluster.pdf
    ln -s $main_dir/${types}/Sample_QC/Cluster*.png Images/${types}_Cluster.png
    ln -s $main_dir/${types}/Sample_QC/cor*.pdf Images/${types}_Corplot.pdf
    ln -s $main_dir/${types}/Sample_QC/cor*.png Images/${types}_Corplot.png
    ln -s $main_dir/${types}/Sample_QC/PCA.pdf Images/${types}_PCA.pdf
    ln -s $main_dir/${types}/Sample_QC/PCA.png Images/${types}_PCA.png
done
cd ..


# ##Differential expression files

# cd DiffExp
# mkdir Images

# for comparison in $comparison_groups
# do
#     for types in bg_lnc bg_m
#     do
# 	# ln -s $main_dir/$types/*.csv .
# 	# ln -s $main_dir/$types/*.xls .
# 	ln -s $main_dir/$types/Diff_exp/Gene_Heatmap.${comparison}*png Images/${types}_${comparison}_Gene_heatmap.png
# 	ln -s $main_dir/$types/Diff_exp/Gene_Heatmap.${comparison}*pdf Images/${types}_${comparison}_Gene_heatmap.pdf
# 	ln -s $main_dir/$types/Diff_exp/transcript_Heatmap.${comparison}*png Images/${types}_${comparison}_transcript_heatmap.png
# 	ln -s $main_dir/$types/Diff_exp/transcript_Heatmap.${comparison}*pdf Images/${types}_${comparison}_transcript_heatmap.pdf
# 	ln -s $main_dir/$types/Diff_exp/Gene_Volcanoplot.${comparison}*png Images/${types}_${comparison}_Gene_volcanoPlot.png
# 	ln -s $main_dir/$types/Diff_exp/Gene_Volcanoplot.${comparison}*pdf Images/${types}_${comparison}_Gene_volcanoPlot.pdf
# 	ln -s $main_dir/$types/Diff_exp/transcript_Volcanoplot.${comparison}*png Images/${types}_${comparison}_transcript_volcanoPlot.png
# 	ln -s $main_dir/$types/Diff_exp/transcript_Volcanoplot.${comparison}*pdf Images/${types}_${comparison}_transcript_volcanoPlot.pdf
# 	mkdir ${types}_${comparison}
# 	cd ${types}_${comparison}
# 	ln -s $main_dir/$types/Diff_exp/*${comparison}.xls .
# 	ln -s $main_dir/$types/Diff_exp/*${comparison}.fa .
# 	cd ..
#     done
# done

# cd ..	

# #GO analysis files
# cd GO
#    ln -s $main_dir/GO/* .
# cd ..

# #Pathway analysis files
# cd Pathway
#    ln -s $main_dir/PATHWAY/* .
# cd ..

# cd $main_dir

# cp -r $main_dir/lncRNA-mRNA/results/* $main_dir/Report/Result/lncRNA-mRNA

# mkdir STAT
# mv algSummary.csv raw_readsSummary.csv trim_readsSummary.csv QC_Stat.txt STAT
