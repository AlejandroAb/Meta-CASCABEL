# Meta-CASCABEL

Snakemake pipeline for assembly and binning of metagenomics reads.

**Current version:** 3.0

The pipeline creates different output files which allow the user to explore the data and results in a simple way, as well as facilitate downstream analysis based on the generated output files.

* Different quality control steps on the reads.
* Taxonomy assessment at different levels
* Read trimming and filtering
* Assembly
* Gene calling (bins, contigs or scaffolds)
* Binning 
* Bin evaluation

## Quick start


**Required input files**

The pipeline is designed to analyze one or more metagenomes.
For each single metagenome you should supply the paired end raw reads:

Forward raw reads (fastq or fastq.gz)
Reverse raw reads (fastq or fastq.gz)

In order to only perform the binning, you can also supply a fasta file containing your assemble. In such case, you also need to supply the raw data.

**Download or clone the repository**

> git clone https://github.com/AlejandroAb/Meta-CASCABEL.git

**Initialize directory structure**

In case you want to process more than one metagenome you need to initialize all the target samples. For this, you need to execute the init_sample.sh script.

The script needs 4 arguments, in the exact same order as follows:

* The name of the project
* The name of the sample
* Absolute path to forward raw reads
* Absolute path to reverse raw reads

Example. Initialize the directory to process two metagenomes (sampleA and sampleB) within a *project* called "test_metagenomes"

```sh
Scripts/init_sample.sh test_metagenomes  sampleA  /path/to/rawdata/sampleA_fw.fastq  /path/to/rawdata/sampleA_rv.fastq
Scripts/init_sample.sh test_metagenomes  sampleB  /path/to/rawdata/sampleB_fw.fastq  /path/to/rawdata/sampleB_rv.fastq
```

Please notice how the project is the same for both commands, while the sample and the reference to the raw data changes.

This steps needs to be executed as many times as the number of samples/metagenomes that the user want to analyze.
At the end you should have a directory structure similar to the following:

```
<PROJECT_NAME>
├── samples
    └── <SAMPLE_NAME_A>
    │   └── rawdata
    │       ├── fw.fastq -> /path/to/raw/data/fw_reads.fq
    │       └── rv.fastq -> /path/ro/raw/data/rv_reads.fq
    └── <SAMPLE_NAME_B>
        └── rawdata
            ├── fw.fastq -> /path/to/raw/data/fw_reads.fq
            └── rv.fastq -> /path/ro/raw/data/rv_reads.fq
```

**Edit configuration file**

<ins>Project name</ins>

```yaml
#------------------------------------------------------------------------------#
#                             Project Name                                     #
#------------------------------------------------------------------------------#
# The name of the project for which the pipeline will be executed. This should#
# be the same name used as the first parameter on init_sample.sh script        #
#------------------------------------------------------------------------------#
PROJECT: "test_metagenomes"
```

If you use the init_script_new.sh for multiple samples, make sure to use the same project name here.

<ins>Samples</ins>

```yaml
#------------------------------------------------------------------------------#
#                               SAMPLES                                        #
#------------------------------------------------------------------------------#
# SAMPLES/Libraries you will like to include on the analysis                   #
# Same sample names used  with init_sample.sh script                           #
# Include all the names between quotes, and comma separated                    #
#------------------------------------------------------------------------------#
SAMPLES: ["sampleA", "sampleB"]
```

In the same way, if the init_script_new.sh for multiple samples/metagenomes was used, make sure to enter
the same sample names, quoted and comma-separated.

Go through the rest of the configuration file and choose your options. 

**Run the pipeline**

*dry run*

 snakemake  --configfile config.yaml  -np

*Run*

 snakemake  --configfile config.yaml  


**Output files structure**

```
<PROJECT>
├── runs
│   └── <RUN>
│       └── <SAMPLE>_data
│           ├── taxonomy  #Output from taxonomy profiling tool
│           │   └── <TAXONOMY_PROFILING>.taxonomy.report
│           ├── trimmed
│           │   ├── qc  #FastQC result for trimmed reads
│           │   ├── read1_paired.fq  #Trimmed reads
│           │   ├── read1_singles.fq
│           │   ├── read2_paired.fq
│           │   └── read2_singles.fq
│           ├── assembly_<ASSEMBLER> 
│           │   ├── contigs.fasta   # Assembly - contigs
│           │   ├── scaffolds.fasta # Assembly - scaffolds (if available)
│           │   └── quast  # Assembly statistics
│           ├── bwa-mem  #Assembly mapping against raw reads
│           │   ├── <ANALYSIS>_<ASSEMBLER>_depth.txt  # depth coverage
│           │   ├── <ANALYSIS>_<ASSEMBLER>_mapped_against_cross-assembly_sorted.bam # bam file
│           │   └── <ANALYSIS>_<ASSEMBLER>_mapped_against_cross-assembly_sorted.flagstat #stats
│           ├── binning #The location for the bins vary per method 
│           │   ├── abundance.<method>.tsv  #Information about the bin abundance per method
│           │   ├── binsanity
│           │   │   └── <ANALYSIS>_<ASSEMBLER>
│           │   │       └── BinSanity-Final-bins  #BinSanity bins folder
│           │   ├── concoct
│           │   │   └── <ANALYSIS>_<ASSEMBLER>    #Concoct bins
│           │   ├── das
│           │   │   └── <ANALYSIS>_<ASSEMBLER>
│           │   │       └──DasOut_DASTool_bins    #DASTool bins
│           │   ├── maxbin
│           │   │   └── <ANALYSIS>_<ASSEMBLER>    #MaxBin bins
│           │   ├── metabat2
│           │   │   └── <ANALYSIS>_<ASSEMBLER>    #Metabat bins
│           │   ├── checkM_<bin_method>
│           │   │   └── summary.txt
│           │   ├── gtdbtk_<bin_method>
│           │   │   ├── gtdbtk.ar122.summary.tsv -> classify/gtdbtk.ar122.summary.tsv
│           │   │   └── gtdbtk.bac120.summary.tsv -> classify/gtdbtk.bac120.summary.tsv
│           │   ├── FinalBins
│           │   │   ├── contig_coverage.txt
│           │   │   ├── new_names.txt
│           │   │   ├── NIOZ114-1.fna
│           │   │   ├── NIOZ114-2.fna
│           │   │   └── NIOZ114-3.fna
│           │   └── FinalBins.summary.tsv
│           └── unbinned
│               ├── unbinned_contigs_list.txt # List of unbinned contigs
│               └── unbinned.fasta # fasta file with unbinned contigs
└── samples
```


**Dependencies** 

* Assembly (you don't need to have all the tools installed, only the one for your target analysis)
  * Spades
  * Megahit
  * IDBA
* Mapping back reads to the assembly
  * BWA
* Binning (you don't need to have all the tools installed, only the one for your target analysis)
  * Maxbin
  * Metabat2
  * CONCOCT
  * Bin Sanity
  * DAS Tool
* Binning evaluation
  * CheckM
  * GTDB-Tk

