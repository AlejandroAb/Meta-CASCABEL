from sys import stdin
import os
contig_N50="0"
contig_L50="0"
contig_n = "0"
scaffold_N50="0"
scaffold_L50="0"
scaffold_n = "0"

def writeOutput(message):
    with open(snakemake.output[0], "a") as tmplog:
            tmplog.write(message)
            tmplog.close()

#read contigs/report.txt and skip first two lines
report_contigs="##Contigs Quast report\n"
with open(snakemake.input[1]) as contigQ:
    i=0;
    for line in contigQ:
        i=i+1
        if i>2 :
            report_contigs += line
        if "N50" in line :
            contig_N50 = line[10:].strip()
        if "L50" in line :
            contig_L50 = line[10:].strip()
        if "# contigs"  in line and i > 10:
            contig_n = line[10:].strip()
    contigQ.close()
#read contigs/report.txt and skip first two lines
report_scaffolds="##Scaffolds Quast report\n"
with open(snakemake.input[0]) as scaffoldQ:
    i=0;
    for line in scaffoldQ:
        i=i+1
        if i>2 :
            report_scaffolds += line
        if "N50" in line :
            scaffold_N50 = line[10:].strip()
        if "L50" in line :
            scaffold_L50 = line[10:].strip()
        if "# contigs"  in line and i > 10:
            scaffold_n = line[10:].strip()
    scaffoldQ.close()



menu = "\033[91m This step validates the assembly \033[0m\n"
menu+="\033[93m Assembler: \033[0m "+ "\033[92m "+str(snakemake.config["ASSEMBLER"])+" \033[0m\n"
menu+="\033[93m # Contigs: \033[0m " + "\033[92m "+contig_n+" \033[0m\n"
menu+="\033[93m Contigs N50: \033[0m " + "\033[92m "+contig_N50+" \033[0m\n"
menu+="\033[93m Contigs L50: \033[0m " + "\033[92m "+contig_L50+" \033[0m\n"
menu+="\033[93m # Scaffolds / Scaffolds broken: \033[0m " + "\033[92m "+scaffold_n+" \033[0m\n"
menu+="\033[93m Scaffolds N50: \033[0m " + "\033[92m "+scaffold_N50+" \033[0m\n"
menu+="\033[93m Scaffolds L50: \033[0m " + "\033[92m "+scaffold_L50+" \033[0m\n"
menu+="\033[91m Please select one option \033[0m\n"
menu+="\033[93m 1. Continue with the work-flow\033[0m\n"
menu+="\033[93m 2. See complete contig report\033[0m\n"
menu+="\033[93m 3. See complete scaffold report\033[0m\n"
menu+="\033[93m 4. Interrupt workflow and delete assembly files \033[0m\n"
menu+="\033[93m 5. Interrupt workflow \033[0m\n"

if snakemake.config["interactive"] == "F":
    print("\033[93m" +"Interactive mode off \033[0m")
    print("\033[93m" +"We suggest to review the full Quast report at: "+ snakemake.input[0]+ " for scaffolds\033[0m")
    print("\033[93m" +"And: "+ snakemake.input[0]+ " for the contigs report\033[0m")
    print("\033[93m" +"Following the contig report:\033[0m")
    print(report_contigs)
    print("\033[93m" +"Following the scaffold report:\033[0m")
    print(report_contigs)
    with open(snakemake.output[0], "w") as tmplog:
        tmplog.write("Interactive mode. Quast validation skipped")
        tmplog.close()
    exit(0)


user_input="0"
while (user_input != "1" and user_input !=  "4" and user_input != "5"):
    if user_input == "2":
        print(report_contigs)
    if user_input == "3":
        print(report_scaffolds)
    print(menu)
    print("\033[92m Enter your option: \033[0m")
    user_input = stdin.readline()
    user_input = user_input[:-1]
    if user_input == "1":
        writeOutput("Option selected: 1. Continue with the work-flow" )
    elif user_input == "4":
        writeOutput("Option selected: 4. Interrupt workflow and delete assembly files\n" )
        print("Cleaning files...")
        shutil.rmtree(snakemake.params[0])
        writeOutput("Files cleanned\n" )
        print("Aborting workflow...")
        exit(1)
    elif user_input == "5":
        writeOutput("Option selected: 5. Interrupt workflow \n" )
        print("Aborting workflow...")
        exit(1)
