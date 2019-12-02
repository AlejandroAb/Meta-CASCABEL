import os
from sys import stdin
with open(snakemake.input[1]) as splitlog:
    median=0
    for line in splitlog:
        if "Median sequence length:" in line:
            try:
                median = int(float((line[line.find(":")+1:])))
                shorts=median-10
                longs=median+10
                splitlog.close()
                break
            except ValueError:
                print("Error trying to cast: "+ line[line.find(":")+1:])
    if median > 0:
        print("\033[91m This step will remove short and long reads \033[0m")
        print("\033[93m Median sequence length: " + str(median) + " \033[0m")
        print("\033[93m More information at: " + snakemake.input[1] + " \033[0m")
        print("\033[93m And: "+ snakemake.params[0] + "/histograms.txt \033[0m")
        print("\033[93m Please enter the option which fits better for your data: \033[0m")
        print("\033[93m 1. Use values from the configuration file: length < "+str(snakemake.config["rm_reads"]["shorts"])+" and length > "+str(snakemake.config["rm_reads"]["longs"])+ "\033[0m")
        print("\033[93m 2. Use values from median +-10: length < " + str(median-10) + " and length > "+ str(median+10)  +" \033[0m")
        print("\033[93m 3. Specify new values! \033[0m")
        print("\033[93m 4. Interrupt workflow \033[0m")
        user_input="0"
        while (user_input != "1" and user_input !=  "2" and user_input != "3" and user_input != "4"):
            print("\033[92m Enter your option: \033[0m")
            user_input = stdin.readline() #READS A LINE
            user_input = user_input[:-1]
        if user_input == "1":
            shorts = snakemake.config["rm_reads"]["shorts"]
            longs = snakemake.config["rm_reads"]["longs"]
        elif user_input == "3":
            ss=-1
            while ss == -1:
                print("\033[92m Please enter the shortest length allowed: \033[0m")
                ui = stdin.readline() #READS A LINE
                ui = ui[:-1]
                try:
                    ss = int(ui)
                    shorts = ss
                except ValueError:
                    print ("Please enter a valid number")
                    ss = -1
            ll=-1
            while ll == -1:
                print("\033[92m Please enter the longest length allowed: \033[0m")
                ui = stdin.readline() #READS A LINE
                ui = ui[:-1]
                try:
                    ll = int(ui)
                    longs = ll
                except ValueError:
                    print ("Please enter a valid number")
                    ll = -1
        elif user_input == "4":
            print("Aborting workflow...")
            exit(1)
        with open(snakemake.output[1], "a") as tmplog:
            tmplog.write(snakemake.input[0] + "\t" + str(shorts) + "\t" + str(longs) + "\n")
            tmplog.close()
        #print("awk '!/^>/ { next } { getline seq } length(seq) > " + str(shorts) + " && length(seq) < " + str(longs) + " { print $0 \"\\n\" seq }' " + snakemake.input[0] + " > "+ snakemake.output[0])
        #os.system("awk '!/^>/ {{ next }} {{ getline seq }} length(seq) >= {config[rm_reads][shorts]} && length(seq) <= {config[rm_reads][longs]} {{ print $0 \"\\n\" seq }}' " + input[0] + " > {output}")
        os.system("awk '!/^>/ { next } { getline seq } length(seq) > " + str(shorts) + " && length(seq) < " + str(longs) + " { print $0 \"\\n\" seq }' " + snakemake.input[0] + " > " + snakemake.output[0])
