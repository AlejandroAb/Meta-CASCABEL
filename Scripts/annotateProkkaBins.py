import os
import sys
import subprocess
#from sys import stdin
#ext = sys.argv[1]
output_dir = snakemake.params["output_dir"]
file_extension = snakemake.params["file_ext"]
print("\033[93mProcessing bin files from directory: \033[0m  \033[92m "+ output_dir+"\033[0m \033[93m with extension:\033[0m \033[92m"+file_extension + " \033[0m")
count = 0
#compute the number of bins to process
for file in os.listdir(output_dir):
    if file.endswith(file_extension):
        count +=1
i = 0
extra_params = []
if len(str(snakemake.config["prokka"]["extra_params"]))>2:
    for param in str(snakemake.config["prokka"]["extra_params"]).split(" "):
        extra_params.append(param.rstrip())

#output_dir+file
with open(snakemake.output[0], "a") as outfile:
    for file in os.listdir(output_dir):
        if file.endswith(file_extension):
            i = i+1
            print("\033[93mAnnotating\033[0m \033[92m" + file  + "\033[0m \033[93m file \033[0m \033[92m" + str(i) + "/" + str(count) + "\033[0m")
            splittedName = file.split(".") #the name is bin.##.extension
            number = splittedName[1]
            fullFIle =[]
            fullFIle.append(output_dir+file)
            prokka = ['prokka',
            "--prefix", "bacteria.bin." +str(number), #Filename output prefix
            "--locustag", "b.bin." +str(number), #Locus tag prefix
            "--species", "Metagenome",
            "--strain", "bacteria.strain." +str(number),
            "--outdir", output_dir+ "prokka_bacteria",
            "--cpus",str(snakemake.config["prokka"]["cpus"]),
            "--addgenes", #Add 'gene' features for each 'CDS' feature
            "--force", #Force overwriting existing output folder
            "--metagenome",# Improve gene predictions for highly fragmented genomes
            "--quiet"  #No screen output
            ]+extra_params+fullFIle
            try:
                status = subprocess.check_call(prokka)
                outfile.write(file + "\tProcessed: OK\n" + str(prokka) + "\n")
            except CalledProcessError:
                print("ERROR " + str(status))
                outfile.write("ERROR " + str(status)+"\n")

            prokka_archaea =['prokka',
            "--prefix", "archaea.bin." +str(number), #Filename output prefix
            "--locustag", "a.bin." +str(number), #Locus tag prefix
            "--species", "Metagenome",
            "--strain", "archaea.strain." +str(number),
            "--outdir", output_dir+ "prokka_archaea",
            "--cpus",str(snakemake.config["prokka"]["cpus"]),
            "--addgenes", #Add 'gene' features for each 'CDS' feature
            "--force", #Force overwriting existing output folder
            "--metagenome",# Improve gene predictions for highly fragmented genomes
            "--kingdom","Archaea", # Annotation mode
            "--quiet"  #No screen output
            ]+extra_params+fullFIle
        #    try:
        #        status = subprocess.check_call(prokka_archaea)
        #        outfile.write(file + "\tProcessed: OK\n" + str(prokka_archaea) + "\n")
        #    except CalledProcessError:
        #        print("ERROR " + str(status))
        #        outfile.write("ERROR " + str(status)+"\n")

    outfile.close()
