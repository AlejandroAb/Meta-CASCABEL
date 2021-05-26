#!/bin/bash
#needs 4 arguments

if  [ $# -lt 3 ]; then
    echo -e "This Script needs 4 arguments\n"
    echo -e "Path to bin folder. e.g. {PROJECT}/runs/{run}/{sample}_data/binning/metabat2/CONTIGS_spades/"
    echo -e "Bins extension. e.g. fa"
    echo -e "PATH to coverage file: e.g. {PROJECT}/runs/{run}/{sample}_data/bwa-mem/CONTIGS_spades_depth.txt"
    echo -e "Output file. e.g.  {PROJECT}/runs/{run}/{sample}_data/binning/abundance.metabat.tsv"
   
   exit 1
fi


for file in $1/*.$2; do
  bin=$(echo $file | awk -F"/" '{print $NF}' | sed 's/\.fa//')
  method=$(echo $file | awk -F"/" '{print $6}');
  grep ">" $file | cut -f2 -d">" |  cut -f1 -d" " | grep -F -w -f - $3 | awk -v bin="$bin" -v m="$method" '{totSize+=$2; totCvg+=$3;}END{print m"\t"bin"\t"NR"\t"totSize"\t"totCvg/NR}' | sort -gk4 >> $4.tmp;
done

echo -e "method\tbin\tnum_contigs\ttot_length\tavg_coverage" | cat -  $4.tmp > $4
rm  $4.tmp
exit 0

