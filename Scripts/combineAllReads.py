#"{PROJECT}/samples/{run}/seqs_fw_rev_filtered.fasta"
#expand("{PROJECT}/samples/{sample}/data/{run}/seqs_fw_rev_filtered.fasta",PROJECT=config["PROJECT"],sample=config["SAMPLES"], run=run)
#This rule was simply executed with the final cat command, but we need the information
#for printing the title in the final report, therefore we create this step, which writedown
#the sample/library names on one file which is readed further on the report_all.py script 
#
import sys
import os
#run = sys.argv[0]
samplesout = snakemake.wildcards.PROJECT + "/runs/" + snakemake.wildcards.run + "/samples.log"
samplesout2 = snakemake.wildcards.PROJECT + "/runs/" + snakemake.wildcards.run + "/cat_samples.log"
allF = snakemake.input.allFiltered
#inputs = allF.split(",")
samples = ""
sampList = ""
print(snakemake.input.allFiltered)
#print(inputs)
for inp in allF:
    sampList += inp + " "
    dirs = inp.split("/")
    if len(dirs) > 2 :
        samples+=dirs[2] + " "
with open(samplesout, "w") as sampfile:
    sampfile.write(samples)
    sampfile.close()
with open(samplesout2, "w") as catsampfile:
    catsampfile.write("cat " + sampList + " > " + snakemake.output[0])
    catsampfile.close()
os.system("cat " + sampList + " > " + snakemake.output[0])
