#!/bin/bash
js() {
    jq ".$1" $config | sed 's/"//g'
}

config="analysis-profile.json"
main_dir="$(js main_dir)"
report_template="$(js report_template)"

echo "report_template=$report_template"

cp -Rf $report_template/* $main_dir/Report
cp $main_dir/STAT/* $main_dir/Report/Tables
rm $main_dir/Report/Tables/QC_Stat.txt

cd $main_dir/Report
python prepare-report.py -t report-template.html -c report-config.json -m report.md > report.html
rm prepare-report.py report-template.html report.md report-config.json
cd $main_dir

# mkdir LOG
# mv *.err LOG
# mv *.out LOG

# curdir=$main_dir
# reportname=$(basename $curdir)
# echo "$reportname"
# zip -r $reportname-report.zip Report
# md5sum $reportname-report.zip > $reportname-report.zip.md5
