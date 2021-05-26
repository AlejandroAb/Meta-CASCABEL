#!/bin/bash
#needs 5 arguments
#$1 Bin directory;
#$2 Bin extension
#$3 SampleName
#$4 Outdir
#$5 Outfile (old to new name)

if  [ $# -lt 5 ]; then
    echo -e "This script needs 5 arguments\n"
   exit 1
fi

bin_num=0;

for file in $1/*.$2; do
#  ((bin_num=bin_num+1));
  bin_num=$((bin_num+1))
  new_name=$3-${bin_num}.fna
  new_name_short=$3-${bin_num}
  old_name=$(echo $file | awk -F"/" '{print $NF}')
  method=$(echo $file | awk -F"/" '{print $6}');
  #cat $file | awk -v  dir="$1" '/>/{sub(">","&"FILENAME"|");sub(/\.fna/,x);sub(,x); >  $4/${new_name}.fna
 # cat $file | sed "s/>/>$new_name_short|/" >  $4/${new_name}
  cat $file | awk -v nm=${new_name_short}   '{if($0 ~ ">"){i=i+1;gsub(">",">"nm".C"i" ",$0);print $0}else{print $0}}' | \
  awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' > $4/${new_name}
  #echo -e  ${method}"\t"${old_name}"\t"${new_name} >> $5
  printf  ${method}"\t"${old_name}"\t"${new_name}"\n" >> $5
done

exit 0

