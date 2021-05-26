
#!/bin/bash
#needs 3 arguments
#$1 Bin directory;
#$2 Bin extension
#$3 Outfile 
#Computes average GC for all the files at $1/*.$2

if  [ $# -lt 3 ]; then
    echo -e "This script needs 3 arguments\n"
   exit 1
fi

for file in $1/*.$2; do
  ext=$2
  #bin_name=$(echo $file | awk -F"/" '{print $NF}' | sed "s/${ext}//\1")
  bin_name=$(echo $file | awk -F"/" '{print $NF}' | awk -F"." '{name=$1; for(i=2;i<NF;i++){name=name"."$i} print name}')
  method=$(echo $file | awk -F"/" '{print $6}');
  perl Scripts/length+GC.pl $file |  awk -F"\t" -v bn=${bin_name} -v m=${method} '{avg_gc=avg_gc+$2}END{print m"\t"bn"\t"avg_gc*100/NR}' >> $3
done

exit 0

