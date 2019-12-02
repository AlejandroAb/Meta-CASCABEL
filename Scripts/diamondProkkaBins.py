import os
import sys
import subprocess
#from sys import stdin
#ext = sys.argv[1]
#https://stackoverflow.com/questions/273192/how-can-i-create-a-directory-if-it-does-not-exist
def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

output_dir = snakemake.params["output_dir"]
prokka_bacteria = output_dir + "prokka_bacteria/"
prokka_archaea = output_dir + "prokka_archaea/"
file_extension = snakemake.params["file_ext"]
output_dir_diamond = snakemake.params["output_dir"] + "diamond/"
ensure_dir(output_dir_diamond)
#print("\033[93mPredicting bin genes from directory: \033[0m  \033[92m "+ prokka_bacteria+"\033[0m \033[93m with extension:\033[0m \033[92m"+file_extension + " \033[0m")
print("\033[93m**** Diamond annotation ****\033[0m")
count = 0
#compute the number of bins to process
for file in os.listdir(prokka_bacteria):
    if file.endswith(file_extension):
        count +=1
i = 0
extra_params = []
out_fmt=[]
if len(str(snakemake.config["diamond"]["extra_params"]))>2:
    for param in str(snakemake.config["diamond"]["extra_params"]).split(" "):
        extra_params.append(param.rstrip())
for opt in str(snakemake.config["diamond"]["outfmt"]).split(" "):
    out_fmt.append(opt.rstrip())

output_dir+file
with open(snakemake.output[0], "a") as outfile:
    for file in os.listdir(prokka_bacteria):
        if file.endswith(file_extension):
            splittedName = file.split(".") #the name is bacteria.bin.##.extension
            number = splittedName[2]
            if not os.path.isfile(output_dir+"diamond/"+"bin."+number+".out"):
                i = i+1
                print("\033[93mData Base: \033[0m \033[92m" + str(snakemake.config["diamond"]["db"])  + "\033[0m")
                print("\033[93mStarting annotation for genes in file: \033[0m \033[92m" + file  + "\033[0m \033[93m file \033[0m \033[92m" + str(i) + "/" + str(count) + "\033[0m")
                diamond= ['diamond',
                "blastp",
                "--db",str(snakemake.config["diamond"]["db"]),
                "--query",prokka_bacteria+file,
                "--out",output_dir+"diamond/"+"bin."+number+".out",
                "--threads",str(snakemake.config["diamond"]["threads"]),
                "--taxonmap",str(snakemake.config["diamond"]["taxonmap"]),
                "--outfmt"#,str(snakemake.config["diamond"]["outfmt"])
                ]+out_fmt+extra_params
                try:
                    status = subprocess.check_call(diamond)
                    outfile.write(file + "\tProcessed: OK\n" + str(diamond) + "\n")
                except CalledProcessError:
                    print("ERROR " + str(status))
                    outfile.write("ERROR " + str(status)+"\n")
            else:
                i =i+1
                print("\033[93mAnnotation already exists for file: \033[0m \033[92m" + file  + "\033[0m \033[93m file \033[0m \033[92m" + str(i) + "/" + str(count) + "\033[0m")
                outfile.write("File: "+ output_dir+"diamond/"+"bin."+number+".out " + "already exists...");

#Prokka uses prodigal for gene calling, and prodigal doent differentiate between bacterial and archaea for gene calling, so at the end, with prokka
#running twice, one for bacterial and one for archaea, we end up with the same group of predicted genes, then prokka performs different annotations
#according to domain, so there is no need to run diamond on both, prokka_bacteria and prokka_archaea or yes? the logic option to do that is that
#prodigal have differente prediction modes  Anonymous (for metagenomes!) NOrmal and Training and only in case that prokka has specialized trainned data
#to call bacterial and archaea  http://redmine.nioz.nl/issues/106
#    for file in os.listdir(prokka_archaea):
#        if file.endswith(file_extension):
#            i = i+1
#            print("\033[93mDiamond annotation for genes in file: \033[0m \033[92m" + file  + "\033[0m \033[93m file \033[0m \033[92m" + str(i) + "/" + str(count) + "\033[0m")
#            splittedName = file.split(".") #the name is bacteria.bin.##.extension
#            number = splittedName[2]
#            diamond= ['diamond',
#            "blastp",
#            "--db",str(snakemake.config["diamond"]["db"]),
#            "--query",output_dir+file,
#            "--out",prokka_bacteria+"diamond/"+"bacteria.bin."+number+".out"
#            "--threads",str(snakemake.config["diamond"]["threads"]),
#            "--taxonmap",str(snakemake.config["diamond"]["taxonmap"])
#            "--outfmt",str(snakemake.config["diamond"]["outfmt"])
#            ]+extra_params
#            try:
#                status = subprocess.check_call(prokka)
#                outfile.write(file + "\tProcessed: OK\n" + str(prokka) + "\n")
#            except CalledProcessError:
#                print("ERROR " + str(status))
#                outfile.write("ERROR " + str(status)+"\n")



    outfile.close()
