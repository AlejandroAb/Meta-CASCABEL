#!/bin/bash
#needs 1 arguments
#$1 Path to summary files. e.g.  {project}/runs/spades_assembly/{sample}_data/binning/gtdbtk_*/summary.txt;

if  [ $# -lt 1 ]; then
    echo -e "This program needs 1 arguments\n"
    echo -e "Path to GTDB Tk summary files in the form: {project}/runs/spades_assembly/{sample}_data/binning/gtdbtk_*/summary.txt"
   exit 1
fi


for file in $1; do
  method=$(echo $file | cut -f6 -d"/" | cut -f2 -d"_"); 
  path=$(echo $file | cut -f1-6 -d"/"); 
  path_out=$(echo $file | cut -f1-5 -d"/");
  echo "Summarizing results: "$method;
  cat ${path}/gtdbtk.*.summary.tsv  | awk -v m="$method" -F"\t" -v idx_aa=16  -v idx_red=18  -v idx_w=19 'BEGIN{OFS="\t";} {if(NR==1){for(i=1;i<=NF;i++){if($i~"_percent"){idx_aa=i}; if($i=="red_value"){idx_red=i}; if($i~"warnings"){idx_w=i}; }} else if($0 !~ "user_genome"){print m,$1,$2,$idx_aa,$idx_red,$idx_w}}' >> ${path_out}/summary_gtdb.tmp;
  # cat ${path}/gtdbtk.*.summary.tsv | grep -v "^user_genome" | awk ...
done

echo -e "method\tuser_genome\tclassification\taa_percent\tred_value\twarnings" | cat -  ${path_out}/summary_gtdb.tmp > ${path_out}/summary_gtdb.tsv
rm  ${path_out}/summary_gtdb.tmp
exit 0

