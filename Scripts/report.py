import subprocess
from snakemake.utils import report

def readBenchmark( benchFile ):
   with open(benchFile) as bfile:
       txt = "**Benchmark info:**\n\n"
       for l in bfile:
         txt += l + "\n\n"
   return txt;

with open(snakemake.input.counts) as counts:
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
            elif "Total" in tmpLine[3]:
                fType = "TOTAL"

            fileName = filePath[-1]
            if fType != "TOTAL":
                countTxt += "* **File**: " + fileName + " **reads**: " + tmpLine[0] + " **type**: :green:`" + fType + "` \n\n"
    #print(countTxt)
fqv = subprocess.run(['fastqc', '--version'], stdout=subprocess.PIPE)
fqVersion = "**" + fqv.stdout.decode('utf-8').strip() + "**"

fq2faVersion = "**NA**"


removeReadsFile = snakemake.wildcards.PROJECT + "/runs/"+snakemake.wildcards.run+"/"+ snakemake.wildcards.sample + "_data/filter.log"
fqBench = readBenchmark(snakemake.wildcards.PROJECT+"/samples/"+snakemake.wildcards.sample+"/qc/fq.benchmark")
pearBench =readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/peared/pear.benchmark")
combineBench =readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/combine_seqs_fw_rev.benchmark")
rmShorLongBench =readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/filter.benchmark")

try:
    pearv = subprocess.run( ["pear -h | grep 'PEAR v'"], stdout=subprocess.PIPE, shell=True)
    pearversion = "**" + pearv.stdout.decode('utf-8').strip() + "**"
except Exception as e:
    pearversion = "Problem reading version"

shorts = str(snakemake.config["rm_reads"]["shorts"])
longs = str(snakemake.config["rm_reads"]["longs"])
with open(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/filter.log") as trimlog:
    for line in trimlog:
        tokens = line.split("\t")
        if len(tokens)>2:
            shorts = tokens[1]
            longs = tokens[2]


#include user description on the report
desc = snakemake.config["description"]
txtDescription = ""
if len(desc) > 0:
    txtDescription = "\n**User description:** "+desc+"\n"
#bcValidationBench =readBenchmark(snakemake.wildcards.PROJECT+"/metadata/bc_validation/"+snakemake.wildcards.sample+"/validation.benchmark")

report("""
Amplicon Analysis Report for Sample: {snakemake.wildcards.sample}
=====================================================================
    .. role:: commd
    .. role:: red
    .. role:: green

This pipeline is designed to run amplicon sequence analysis.

The main idea is to take raw data and perform all the necessary steps, in order to create different output files which allow the user to explore it's data on a simply and meaningful way, as well as facilitate downstream analysis, based on the generated output files.
{txtDescription}
One of the goals of this pipeline is also to encourage the documentation in order to help with the data analysis reproducibility.

Following you can see all the steps that were taken in order to get this pipeline's final results.

Raw Data
---------
The raw data for this libraries can be found at:

:green:`- FW raw reads:` {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.sample}/rawdata/fw.fastq

:green:`- RV raw reads:` {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.sample}/rawdata/rv.fastq

FastQC
-------
FastQC is a program to evaluate the sequence quality.

This step was performed with the following command:

:commd:`fastqc {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.sample}/rawdata/fw.fastq {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.sample}/rawdata/rv.fastq --extract -o {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.sample}/qc/`

:red:`The version used for this run was:` {fqVersion}

You can follow the link in order to see complete FastQC report:

FastQC for sample {snakemake.wildcards.sample}_1: FQ1_
    .. _FQ1: ../qc/fw_fastqc.html

FastQC for sample {snakemake.wildcards.sample}_2: FQ2_
    .. _FQ2: ../qc/rv_fastqc.html

{fqBench}

PEAR
-----
PEAR is a fast and accurate Illumina Paired-End reAd mergeR. This tool is used to extend paired end fragments

This step was performed with the following command:

:commd:`pear -f {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.sample}/rawdata/fw.fastq -r {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.sample}/rawdata/rv.fastq -t {snakemake.config[pear][t]} -v {snakemake.config[pear][v]} -j {snakemake.config[pear][j]} -p {snakemake.config[pear][p]} -o {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/peared/seqs > {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/peared/seqs.assembled.fastq`

:red:`The version used for this run was:` {pearversion}

**The output files for this step are:**

:green:`- Merged reads:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/peared/seqs.assembled.fastq

:green:`- Log file:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/peared/pear.log

{pearBench}

Extract barcodes
----------------
A variety of data formats are possible, depending upon how one utilized sequencing primers, designed primer constructs (e.g., partial barcodes on each end of the read), or processed the data (e.g., barcodes were put into the sequence labels rather than the reads)

On that regard, this step is performed in order to extract the barcodes used for primer sequencing.

This script is designed to format fastq sequence and barcode data so they are compatible with split_libraries_fastq.py

This step was performed with the following command:

:commd:`extract_barcodes.py -f {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/peared/seqs.assembled.fastq -c {snakemake.config[ext_bc][c]}  --bc1_len {snakemake.config[ext_bc][bc1_len]} --bc2_len {snakemake.config[ext_bc][bc2_len]} -o {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/barcodes/`

:red:`The version used for this run was:` {ebVersion}

**The output files for this step are:**

:green:`- Fastq file with barcodes:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/barcodes/barcodes.fastq

:green:`- Fastq file with the reads:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/barcodes/reads.fastq

{barBench}

{correctBCStr}

Split libraries
----------------
This step performs demultiplexing of Fastq sequence data where barcodes and sequences are contained in two separate fastq files (common on Illumina runs)

The command run twise, one for forward reads and one for reverse reads and was performed with the following command:

:commd:`split_libraries_fastq.py -m {snakemake.wildcards.PROJECT}/metadata/sampleList_mergedBarcodes_{snakemake.wildcards.sample}.txt -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/barcodes/reads.fastq -o  {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/splitLibs -b {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/barcodes/{bcFile} -q {snakemake.config[split][q]} --barcode_type {snakemake.config[split][barcode_type]} {snakemake.config[split][extra_params]} && split_libraries_fastq.py -m {snakemake.wildcards.PROJECT}/metadata/sampleList_mergedBarcodes_{snakemake.wildcards.sample}.txt -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/barcodes/reads.fastq -o {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/splitLibsRC -b {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/barcodes/{bcFile} -q {snakemake.config[split][q]} --barcode_type {snakemake.config[split][barcode_type]} {snakemake.config[split][extra_params]} --rev_comp_mapping_barcodes --rev_comp`

:red:`The version used for this run was:` {spVersion}

**The output files for this step are:**

:green:`- FW reads fasta file with new header:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/splitLibs/seqs.fna

:green:`- Text histogram with the length of the fw reads:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/splitLibs/histograms.txt

:green:`- Log file for the fw reads:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/splitLibs/split_library_log.txt

:green:`- RV reads fasta file with new header:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/splitLibsRC/seqs.fna

:green:`- Text histogram with the length of the rv reads:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/splitLibsRC/histograms.txt

:green:`- Log file for the rv reads:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/splitLibsRC/split_library_log.txt

{splitLibsBench}

Combine reads
--------------
The next step consists on the creation of one single file with the forward and reverse reads.

This operation is performed with the cat command:

:commd:`cat {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/splitLibs/seqs.fna {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/splitLibsRC/seqs.fna > {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/seqs_fw_rev_accepted.fna`

**The output file for this step is:**

:green:`- Fasta file with combined reads:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/seqs_fw_rev_accepted.fna

{combineBench}

Remove too long and too short reads
------------------------------------
At this point we already have a single file, but before continue to the OTU classification is necessary to remove short and long reads.

In this case the executed command is:

:commd:`awk '!/^>/ {{ next }} {{ getline seq }} length(seq) > shorts  && length(seq) < longs {{ print $0 \"\\n\" seq }}'  {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/seqs_fw_rev_accepted.fna  >  {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/seqs_fw_rev_filtered.fasta`

**The output file for this step is:**

:green:`- Fasta file with correct sequence length:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/seqs_fw_rev_filtered.fasta

{rmShorLongBench}

Final counts
-------------

{countTxt}


Taxonomic combined report
---------------------------
The taxonomic report performed in combine with the different supplied libraries, could be find at the following link: taxo_report_
    .. _taxo_report: ../../../runs/{snakemake.wildcards.run}/report_all.html


""", snakemake.output[0], metadata="Author: J. Engelmann & A. Abdala ")
