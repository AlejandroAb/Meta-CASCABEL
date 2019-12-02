import os
from sys import stdin
htmlFile="Not Found"
logFile="Not Found"
for file in os.listdir(snakemake.params[0]):
    if file.endswith("."):
        logFile=snakemake.params[0]+file
    elif file.endswith(".html"):
        htmlFile=snakemake.params[0]+file
with open(logFile) as bcvlog:
    for line in bcvlog:
        if not "No errors or warnings found in mapping file" in line:
            print("\033[91m" + "Validation mapping file contains some warnings or errors: " + logFile + "\033[0m")
            print("Please take a look on complete report at: "+ htmlFile)
            print("\033[93m" +"If continue, maybe an error will be thrown during extract_bc rule. Do you want to continue anyway y/n?"+ "\033[0m")
            user_input = stdin.readline() #READS A LINE
            user_input = user_input[:-1]
            if user_input.upper() == "Y" or user_input.upper() == "YES":
                print("\033[92m" +"The flow goes on!"+ "\033[0m")
                with open(snakemake.output[0], "w") as tmplog:
                    tmplog.write("Error on barcode validation mapping, user continue...")
                    tmplog.close()
                break
            else:
                print("Aborting workflow...")
                logfile.close()
                exit(1)
        else:
            with open(snakemake.output[0], "w") as tmplog:
                tmplog.write("Barcode validation log OK")
                tmplog.close()
