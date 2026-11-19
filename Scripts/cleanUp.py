
import os
import glob
import fnmatch
import shutil


def remove_file(file, log):
    """
    1. Checks if file exists
    2. Checks if input is a file or directory
    3. Removes file/directory
    4. Writes removed file to log
    """
    if not os.path.lexists(file):
        return
    if os.path.islink(file) or os.path.isfile(file):
        os.remove(file)
    elif os.path.isdir(file):
        shutil.rmtree(file)
    else:
        return
    log.write(f"Removed: {file}\n")

def remove_files(file_path, log):

    # MEGAHIT
    log.write("Starting MEGAHIT cleanup\n")
    megahit_dir = os.path.join(file_path, "assembly_MEGAHIT")
    for name in ["checkpoints.txt", "done"]:
       remove_file(os.path.join(megahit_dir, name), log)
    intermediate_dir = os.path.join(megahit_dir, "intermediate_contigs")
    for file in glob.glob(os.path.join(intermediate_dir, "*")):
        if ".final." not in os.path.basename(file):
            remove_file(file, log)

    # SPADES
    log.write("Starting SPADES cleanup\n")
    """
    Nog even terugkomen op de K/ situationne
    """

    # SPADES
    log.write("Starting SPADES cleanup\n")

    spades_dir = os.path.join(file_path, "assembly_SPADES")
    spades_keep_list = ["*.benchmark","*.fasta.gz","*.txt","*.log","quast","quast/*","K*/final_contigs.fasta","K*/scaffolds.fasta"]

    for file in glob.glob(os.path.join(spades_dir, "**", "*"), recursive=True):
        relative_path = os.path.relpath(file, spades_dir)
        keep = any(fnmatch.fnmatch(relative_path, pattern) for pattern in spades_keep_list)
        if keep:
            continue
        if os.path.isdir(file):
            if any(fnmatch.fnmatch(os.path.relpath(f, spades_dir), pattern) for f in glob.glob(os.path.join(file, "**", "*"), recursive=True) for pattern in spades_keep_list):
                continue
        remove_file(file, log)


    """
    spades_dir = os.path.join(file_path, "assembly_SPADES")
    spades_keep_list = ["*.benchmark","*.fasta.gz","*.txt","*.log","quast","K*/final.contigs.fasta","K*/scaffolds.fasta"]
    for file in glob.glob(os.path.join(spades_dir, "**", "*"),recursive=True):
        relative_path = os.path.relpath(file, spades_dir)
        if relative_path == "quast" or relative_path.startswith("quast" + os.sep):
          continue
        keep = any(fnmatch.fnmatch(relative_path, pattern)for pattern in spades_keep_list)
        if not keep:
            remove_file(file, log)
    """

    # IDBA
    log.write("Starting IDBA cleanup\n")
    idba_dir = os.path.join(file_path, "assembly_IDBA")
    idba_keep_list = ["*.benchmark","*.fasta.gz","*.txt","log","quast"]
    for file in glob.glob(os.path.join(idba_dir, "*")):
        relative_path = os.path.relpath(file, idba_dir)
        keep = any(fnmatch.fnmatch(relative_path, pattern) for pattern in idba_keep_list)
        if not keep:
            remove_file(file, log)

    # Binning
    binning = snakemake.config["BINNING"]

    # Metabat
    log.write("Starting METABAT cleanup\n")
    metabat_dir = os.path.join(file_path, "binning", "metabat2")
    if binning not in ["METABAT", "DAS"]:
        remove_file(metabat_dir, log)

    # Binsanity
    log.write("Starting BINSANITY cleanup\n")
    binsanity_dir = os.path.join(file_path, "binning", "binsanity")
    if binning == "BINSANITY" or (binning == "DAS" and snakemake.config["das"]["binsanity"]["run"] == "T"):
        for file in glob.glob(os.path.join(binsanity_dir, "*", "Binsanity-records*")):
            remove_file(file, log)
    else:
        remove_file(binsanity_dir, log)

    # concoct
    log.write("Starting CONCOCT cleanup\n")
    concoct_dir = os.path.join(file_path, "binning", "concoct")
    if binning == "CONCOCT" or (binning == "DAS" and snakemake.config["das"]["concoct"]["run"] == "T"):
        concoct_remove = os.path.join(concoct_dir, "*", "*.csv")
        for file in glob.glob(concoct_remove):
            remove_file(file, log)
    else:
        remove_file(concoct_dir, log)

    # MaxBin
    log.write("Starting MAXBIN cleanup\n")
    maxbin_dir = os.path.join(file_path, "binning", "maxbin")
    if binning == "MAXBIN" or (binning == "DAS" and snakemake.config["das"]["maxbin"]["run"] == "T"):
        maxbin_keep_list = ["binTable.tsv","*.fasta.gz","*.log","bin.summary",]
        for file in glob.glob(os.path.join(maxbin_dir, "*", "*")):
            if not os.path.isfile(file):
                continue
            filename = os.path.basename(file)
            keep = any(fnmatch.fnmatch(filename, pattern) for pattern in maxbin_keep_list)
            if not keep:
                remove_file(file, log)
    else:
        remove_file(maxbin_dir, log)

    # semibin2
    log.write("Starting SEMIBIN cleanup\n")
    semibin_dir = os.path.join(file_path, "binning", "semibin2")
    if binning == "SEMIBIN" or (binning == "DAS" and snakemake.config["das"]["semibin"]["run"] == "T"):
        for file in glob.glob(os.path.join(semibin_dir, "*", "model.pt")):
            remove_file(file, log)
        semibin_remove = os.path.join(semibin_dir, "*", "*.csv")
        for file in glob.glob(semibin_remove):
            remove_file(file, log)
    else:
        remove_file(semibin_dir, log)

def main():
    with open(snakemake.output[0], "w") as log:
        if snakemake.config["remove_intermediate_output"] == "F":
            log.write("No intermediate files were removed\n")
            log.write("Config file value: "+snakemake.config["remove_intermediate_output"]+"\n")
        else:
            dir_path=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/"
            
            remove_files(dir_path, log)
            log.write("Cleanup completed\n")

main()

