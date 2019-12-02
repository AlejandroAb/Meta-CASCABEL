import os
import sys
import subprocess
import re
#from sys import stdin
#ext = sys.argv[1]
output_dir_metabat = snakemake.params["output_dir_metabat"]
output_dir_maxbin = snakemake.params["output_dir_maxbin"]
output_dir_binsanity = snakemake.params["output_dir_binsanity"]
file_extension_metabat = snakemake.params["file_ext_metabat"]
file_extension_maxbin = snakemake.params["file_ext_maxbin"]
file_extension_binsanity = snakemake.params["file_ext_binsanity"]
concoct_clustering = snakemake.input["concoct"]
bin_sanity_low_completion = snakemake.config["binsanity"]["low_completion"]
#print("\033[93mProcessing bin files from directory: \033[0m  \033[92m "+ output_dir_metabat+"\033[0m \033[93m with extension:\033[0m \033[92m"+file_extension + " \033[0m")

#extra_params = []
#if len(str(snakemake.config["prokka"]["extra_params"]))>2:
#    for param in str(snakemake.config["prokka"]["extra_params"]).split(" "):
#        extra_params.append(param.rstrip())

#output_dir+file
#metabat
for file in os.listdir(output_dir_metabat):
    if file.endswith(file_extension_metabat):
        #i = i+1
        #print("\033[93mAnnotating\033[0m \033[92m" + file  + "\033[0m \033[93m file \033[0m \033[92m" + str(i) + "/" + str(count) + "\033[0m")
        splittedName = file.split(".") #the name is bin.##.extension
        number = splittedName[1]
        fullFile=output_dir_metabat+"/"+file
        os.system("cat "+ fullFile + " | grep '^>' | sed 's/>//g' | awk '{print $1\"\\tmetabat."+number+"\"}' >> " +snakemake.output["metabat_out"])
#maxbin
for file in os.listdir(output_dir_maxbin):
    if file.endswith(file_extension_maxbin):
        #i = i+1
        #print("\033[93mAnnotating\033[0m \033[92m" + file  + "\033[0m \033[93m file \033[0m \033[92m" + str(i) + "/" + str(count) + "\033[0m")
        splittedName = file.split(".") #the name is bin.##.extension
        number = splittedName[1]
        fullFile=output_dir_maxbin+"/"+file
        os.system("cat "+ fullFile + " | grep '^>' | sed 's/>//g' | awk '{print $1\"\\tmaxbin."+number+"\"}' >> " +snakemake.output["maxbin_out"])
#BinSanity
for file in os.listdir(output_dir_binsanity):
    if file.endswith(file_extension_binsanity) and (file.startswith("final") or (bin_sanity_low_completion == "T" and file.startswith("low_"))):
        name = "l" if file.startswith("low_") else "f"
        splittedName = re.split('\.|_|-',file)
        if name.startswith("f") and "refined" in file:
            number = splittedName[-2]
            number1 = splittedName[2]
            name+="."+number1+"."+number
        else:
            number = splittedName[-2] #name is low_completion-refined_17.fna or final_Bin-1670.fna  final_Bin-1670-refined_1
            name+="."+number


        fullFile=output_dir_binsanity+"/"+file
        os.system("cat "+ fullFile + " | grep '^>' | sed 's/>//g' | awk '{print $1\"\\tbinsanity."+name+"\"}' >> " +snakemake.output["binsanity_out"])
#CONCOCT
os.system("cat "+ concoct_clustering + " | awk -F\",\" '{print $1\"\tconcoct.\"$2}'   > " +snakemake.output["concoct_out"])
