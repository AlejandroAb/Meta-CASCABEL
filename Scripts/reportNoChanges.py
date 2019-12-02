import subprocess
from snakemake.utils import report
with open(snakemake.input[0]) as counts:
    #n_calls = sum(1 for l in vcf if not l.startswith("#"))
    countTxt="Following you can see the final counts: \n\n"
    fqVersion = ""
    for l in counts:
        tmpLine = l.split(' ')
        #print(tmpLine)
        if len(tmpLine)>1:
            filePath = tmpLine[3].split("/")
            fType = "undefined"
            if "raw" in tmpLine[3]:
                fType = "raw data"
            elif "fw_rev" in tmpLine[3]:
                fType = "accepted clean reads"
            elif "assembled" in tmpLine[3]:
                fType = "Peared file"
            elif "split" in tmpLine[3]:
                fType = "Split library file"

            fileName = filePath[-1]
            countTxt += "* **File**: " + fileName + " **reads**: " + tmpLine[0] + " **type**: :red:`" + fType + "` \n\n"
    #print(countTxt)
    fqv = subprocess.run(['fastqc', '--version'], stdout=subprocess.PIPE)
    fqVersion = "**" + fqv.stdout.decode('utf-8').strip() + "**"
    ebv = subprocess.run(['extract_barcodes.py', '--version'], stdout=subprocess.PIPE)
    ebVersion = ebv.stdout.decode('utf-8')
    ebVersion = "**" + ebVersion[ebVersion.find(":")+1:].strip() + "**"
    spv = subprocess.run(['split_libraries_fastq.py', '--version'], stdout=subprocess.PIPE)
    spVersion = spv.stdout.decode('utf-8')
    spVersion = "**" + spVersion[spVersion.find(":")+1:].strip() + "**"
    correctBCStr = ""
    bcFile="barcodes.fastq"
    if config["bc_missmatch"]:
        correctBCStr = "Correct Barcode\n------------------\n"
        correctBCStr += "Some times the firsts or lasts bases on reads sequences, tend to have some errors which are represented with N's. This error could lead to dismiss a lot of sequences since we dont know the origin sample of those sequences. In order to deal with this issue, this step try to correct those sequences by assigning them a corrected barcode.\n\n"
        correctBCStr += "The maxim number of missmatches this step corrected is **"  + str(snakemake.config["bc_missmatch"]) + "**. For that, the following *in-house* made R script was executed:\n\n"
        correctBCStr += ":commd:`Rscript Scripts/errorCorrectBarcodes.R $PWD "+snakemake.wildcards.PROJECT+"/metadata/sampleList_mergedBarcodes_"+snakemake.wildcards.sample+".txt "+wildcards.PROJECT+"/samples/"+wildcards.sample+"/data/"+wildcards.run+"/barcodes/barcodes.fastq "  + str(config["bc_missmatch"]) + "`\n\n"
        correctBCStr += "The file with the corrected barcodes can be find at: "+snakemake.wildcards.PROJECT+ "/samples/"+snakemake.wildcards.sample+"/data/"+snakemake.wildcards.run+"/barcodes/barcodes.fastq_corrected\n\n"
        bcFile="barcodes.fastq_corrected"
    removeReadsFile = wildcards.PROJECT + "/samples/"+ wildcards.sample + "/data/"+wildcards.run+"/filter.log"
report("""
Amplicon Analysis Report for Sample: {wildcards.sample}
=====================================================================
    .. role:: commd
    .. role:: red

This pipeline is designed to run amplicon sequence analysis.

The main idea is to take raw data and perform all the necessary steps, in order to create different output files which allow the user to explore it's data on a simply and meaningful way, as well as facilitate further analysis, based on the generated output files.

One of the goals of this pipeline is also to encourage the documentation in order to help with the data analysis reproducibility.

Following you can see all the steps that were taken in order to get this pipeline's final results.

Raw Data
---------
The pipeline main (and mandatory) input consists on the raw data, this can be found at:

**If you want to download this data please rigth click the link and choose save as..**

FW raw reads {wildcards.sample}_1: FW_
    .. _FW: ../data/rawdata/fw.fastq
RV raw reads {wildcards.sample}_2: RV_
    .. _RV: ../data/rawdata/rv.fastq

..    .. _RV: {wildcards.PROJECT}/samples/{wildcards.sample}/data/rawdata/rv.fastq

FastQC
-------
FastQC is a program to evaluate the sequence quality.

The version used for this run was: {fqVersion}

This step was performed with the following command:

:commd:`fastqc {wildcards.PROJECT}/samples/{wildcards.sample}/data/rawdata/fw.fastq {wildcards.PROJECT}/samples/{wildcards.sample}/data/rawdata/rv.fastq --extract -o {wildcards.PROJECT}/samples/{wildcards.sample}/qc/`

You can follow the link in order to see complete FastQC report:

FastQC for sample {wildcards.sample}_1: FQ1_
    .. _FQ1: ../qc/fw_fastqc.html

FastQC for sample {wildcards.sample}_2: FQ2_
    .. _FQ2: ../qc/rv_fastqc.html

PEAR
-----
PEAR is a fast and accurate Illumina Paired-End reAd mergeR. This tool is used to extend paired end fragments

The version used for this run was: **PEAR v0.9.8**

This step was performed with the following command:

:commd:`pear -f {wildcards.PROJECT}/samples/{wildcards.sample}/data/rawdata/fw.fastq -r {wildcards.PROJECT}/samples/{wildcards.sample}/data/rawdata/rv.fastq -t {config[pear][t]} -v {config[pear][v]} -j {config[pear][j]} -p {config[pear][p]} -o {wildcards.PROJECT}/samples/{wildcards.sample}/data/{run}/peared/seqs > {wildcards.PROJECT}/samples/{wildcards.sample}/data/{run}/peared/seqs.assembled.fastq`

The full log is available to download: pear.log_
    .. _pear.log: ../data/{run}/peared/pear.log

The assembled reads can be found at: {wildcards.PROJECT}/samples/{wildcards.sample}/data/{run}/peared/seqs.assembled.fastq

Extract barcodes
----------------
A variety of data formats are possible, depending upon how one utilized sequencing primers, designed primer constructs (e.g., partial barcodes on each end of the read), or processed the data (e.g., barcodes were put into the sequence labels rather than the reads)

On that regard, this step is performed in order to extract the barcodes used for primer sequencing.

The version used for this run was: {ebVersion}

This step was performed with the following command:

:commd:`extract_barcodes.py -f {wildcards.PROJECT}/samples/{wildcards.sample}/data/{run}/peared/seqs.assembled.fastq -c {config[ext_bc][c]}  --bc1_len {config[ext_bc][bc1_len]} --bc2_len {config[ext_bc][bc2_len]} -o {wildcards.PROJECT}/samples/{wildcards.sample}/data/{run}/barcodes/`

{correctBCStr}
Split libraries
----------------
Once that we have the barcodes to identify the sample origin for each read, it is needed to demultiplex this reads by spliting the libraries.

The version used for this run was: {spVersion}

This step was performed with the following command:

:commd:`split_libraries_fastq.py -m {wildcards.PROJECT}/metadata/sampleList_mergedBarcodes_{wildcards.sample}.txt {wildcards.PROJECT}/samples/{wildcards.sample}/data/{run}/barcodes/reads.fastq -o  {wildcards.PROJECT}/samples/{wildcards.sample}/data/{run}/splitLibs -b {wildcards.PROJECT}/samples/{wildcards.sample}/data/{run}/barcodes/{bcFile} -q {config[split][q]} --barcode_type {config[split][barcode_type]} {config[split][extra_params]} && split_libraries_fastq.py -m {wildcards.PROJECT}/metadata/sampleList_mergedBarcodes_{wildcards.sample}.txt {wildcards.PROJECT}/samples/{wildcards.sample}/data/{run}/barcodes/reads.fastq -o {wildcards.PROJECT}/samples/{wildcards.sample}/data/{run}/splitLibsRC -b {wildcards.PROJECT}/samples/{wildcards.sample}/data/{run}/barcodes/{bcFile} -q {config[split][q]} --barcode_type {config[split][barcode_type]} {config[split][extra_params]} --rev_comp_mapping_barcodes --rev_comp`

The command run twise, one for forward reads and one for reverse reads.

The log for both runs can be found at: split.fw.log_ and split.rv.log_
    .. _split.fw.log: ../data/{run}/splitLibs/split_library_log.txt
    .. _split.rv.log: ../data/{run}/splitLibsRC/split_library_log.txt

Combine reads
--------------
The next step consists on the creation of one single file with the forward and reverse reads.

This operation is performed with the cat command:

:commd:`cat {wildcards.PROJECT}/samples/{wildcards.sample}/data/{run}/splitLibs/seqs.fna {wildcards.PROJECT}/samples/{wildcards.sample}/data/{run}/splitLibsRC/seqs.fna > {wildcards.PROJECT}/samples/{wildcards.sample}/data/{run}/seqs_fw_rev_accepted.fna`

Remove too long and too short reads
------------------------------------
At this point we already have a single file, but before continue to the OTU classification is necessary to remove short and long reads.

In this case the parameters where:

Final counts
-------------

{countTxt}

More infor for the counts: counts_

The end
""", output[0], metadata="Author: J. Engelmann & A. Abdala ", counts=input[0])
