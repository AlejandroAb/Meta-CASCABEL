#!/bin/bash
#needs 1 arguments
#$1 Path to summary files. e.g.  {project}/runs/spades_assembly/{sample}_data/binning/checkM_*/summary.txt;

if  [ $# -lt 1 ]; then
    echo -e "This program needs 1 arguments\n"
    echo -e "Path to GTDB Tk summary files in the form: {project}/runs/spades_assembly/{sample}_data/binning/checkM_*/summary.txt"
   exit 1
fi


for file in $1; do
  method=$(echo $file | cut -f6 -d"/" | cut -f2 -d"_");
  path_out=$(echo $file | cut -f1-5 -d"/");
  echo $method;
  cat $file | grep -v "^-" | sed 's/ (/_(/'| awk -v m=$method 'BEGIN{OFS="\t"} NR>1{$1=m"\t"$1;print}' >> ${path_out}/summary_checkM.tmp
done

echo -e "method\tBin Id\tMarker_lineage\t#genomes\t#markers\t#marker_sets\t0\t1\t2\t3\t4\t5+\tCompleteness\tContamination\tStrain_heterogeneity" | cat -  ${path_out}/summary_checkM.tmp >  ${path_out}/summary_checkM.tsv
rm  ${path_out}/summary_checkM.tmp
exit 0

