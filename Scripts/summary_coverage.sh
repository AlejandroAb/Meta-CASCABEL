#!/bin/bash
#needs 1 arguments
#$1 Path to summary files. e.g.  {project}/runs/spades_assembly/{sample}_data/binning/gtdbtk_*/summary.txt;

if  [ $# -lt 1 ]; then
    echo -e "This program needs 1 arguments\n"
    echo -e "Path to bin folder: {project}/runs/spades_assembly/{sample}_data/binning/gtdbtk_*/summary.txt"
   exit 1
fi


for file in $1; do
  method=$(echo $file | cut -f6 -d"/" | cut -f2 -d"_"); 
  path=$(echo $file | cut -f1-6 -d"/"); 
  path_out=$(echo $file | cut -f1-5 -d"/");
  echo $method; 
  cat ${path}/gtdbtk.*.summary.tsv | grep -v "^user_genome" | awk -v m=$method -F"\t" 'BEGIN{OFS="\t"}{print m,$1,$2,$16,$18}' >> ${path_out}/summary_gtdb.tmp;
done

echo -e "method\tuser_genome\tclassification\taa_percent\tred_value" | cat -  ${path_out}/summary_gtdb.tmp > ${path_out}/summary_gtdb.tsv
rm  ${path_out}/summary_gtdb.tmp
exit 0

