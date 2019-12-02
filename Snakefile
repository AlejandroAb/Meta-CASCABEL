"""
Metagenomics Workflow for NIOZ MMBL.
Version: 2.1
Author: Julia Engelmann and Alejandro Abdala
Last update: 24/10/2018
"""
run=config["RUN"]
rule all:
    input:
        expand("{PROJECT}/runs/{run}/{sample}_data/report_f.html", PROJECT=config["PROJECT"],sample=config["SAMPLES"], run=run)
#        "{PROJECT}/runs/{run}/{sample}_data/report_f.html"
        #mapped_against_cross-assembly_sorted.flagstat
if len(config["SAMPLES"])==1 and len(config["fw_reads"])>0 and len(config["rv_reads"])>0:
    rule init_structure:
        input:
            fw = config["fw_reads"],
            rv = config["rv_reads"]
        output:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq",
            r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq"
        shell:
            "./init_sample.sh "+config["PROJECT"]+" "+config["SAMPLES"][0]+"  {input.fw} {input.rv}"


if config["FastQC"]["onRawReads"] == "T":
#First we run fastQC over the rawdata
#First we run fastQC over the rawdata
    rule fast_qc:
        input:
            r1="{PROJECT}/samples/{sample}/rawdata/fw.fastq",
            r2="{PROJECT}/samples/{sample}/rawdata/rv.fastq"
        output:
            o1="{PROJECT}/samples/{sample}/qc/fw_fastqc.html",
            o2="{PROJECT}/samples/{sample}/qc/rv_fastqc.html",
            s1="{PROJECT}/samples/{sample}/qc/fw_fastqc/summary.txt",
            s2="{PROJECT}/samples/{sample}/qc/rv_fastqc/summary.txt"
        benchmark:
            "{PROJECT}/samples/{sample}/qc/fq.benchmark"
        shell:
            "fastqc {input.r1} {input.r2} --extract -o {wildcards.PROJECT}/samples/{wildcards.sample}/qc/"
    #validate qc if to many fails on qc report
    rule validateQC:
        input:
            "{PROJECT}/samples/{sample}/qc/fw_fastqc/summary.txt",
            "{PROJECT}/samples/{sample}/qc/rv_fastqc/summary.txt",
            "{PROJECT}/samples/{sample}/qc/fw_fastqc.html",
            "{PROJECT}/samples/{sample}/qc/rv_fastqc.html",
            "{PROJECT}/samples/{sample}/rawdata/fw.fastq",
            "{PROJECT}/samples/{sample}/rawdata/rv.fastq"
        output:
            "{PROJECT}/samples/{sample}/qc/fq_fw_internal_validation.txt",
            "{PROJECT}/samples/{sample}/qc/fq_rv_internal_validation.txt"
        script:
            "Scripts/validateQC.py"
else:
    rule skip_raw_fastqc:
        output:
            fw="{PROJECT}/samples/{sample}/qc/fq_fw_internal_validation.txt",
            rv="{PROJECT}/samples/{sample}/qc/fq_rv_internal_validation.txt"
        shell:
            "echo \"User skip raw reads quality control\" > {output.fw} && echo \"User skip raw reads quality control\" > {output.rv}"
#Run LowQA
#rule lowQuality:
#    input:
#        "{PROJECT}/samples/{sample}/rawdata/fw.fastq",
#        "{PROJECT}/samples/{sample}/rawdata/rv.fastq",
#        "{PROJECT}/samples/{sample}/qc/fq_fw_internal_validation.txt",
#        "{PROJECT}/samples/{sample}/qc/fq_rv_internal_validation.txt"
#    output:
#        ""

# run trimmomatic to remove adapter contamination and trim very low quality parts (ends) of the reads.
# trimmomatic-0.35.jar PE -threads 2 $inFile1 $inFile2 $fol/read1_paired.fq $fol/read1_singles.fq $fol/read2_paired.fq $fol/read2_singles.fq
# ILLUMINACLIP:/usr/local/bioinf/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:12:2:TRUE MAXINFO:40:0.6 MINLEN:40
rule trimmomatic:
    input:
        fw="{PROJECT}/samples/{sample}/rawdata/fw.fastq",
        rv="{PROJECT}/samples/{sample}/rawdata/rv.fastq",
        tmp1="{PROJECT}/samples/{sample}/qc/fq_fw_internal_validation.txt",
        tmp2="{PROJECT}/samples/{sample}/qc/fq_rv_internal_validation.txt"
    output:
        read1_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
        read1_single="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_singles.fq",
        read2_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq",
        read2_single="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_singles.fq"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/trimmed/trimmomatic.benchmark"
    shell:
        "java -jar /opt/biolinux/Trinity-v2.6.6/trinity-plugins/Trimmomatic-0.36/trimmomatic-0.36.jar {config[trimm][mode]} -threads {config[trimm][threads]} {input.fw} {input.rv} "
        "{output.read1_paired} {output.read1_single} {output.read2_paired} {output.read2_single} "
        "{config[trimm][clip][type]}:{config[trimm][clip][adapter]}:{config[trimm][clip][seed]}:{config[trimm][clip][palindrome_ct]}:"
        "{config[trimm][clip][simple_ct]}:{config[trimm][clip][minAdpLength]}:{config[trimm][clip][keepBoth]} "
        "{config[trimm][maxinfo][type]}:{config[trimm][maxinfo][targetLength]}:{config[trimm][maxinfo][strictness]} "
        "{config[trimm][minlen][type]}:{config[trimm][minlen][len]}"
if config["FastQC"]["onTrimmedReads"] == "T":
    rule qc_trimmed_reads:
        input:
            r1="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
            r2="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq"
        output:
            o1="{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/read1_paired.fq_fastqc.html",
            o2="{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/read2_paired.fq_fastqc.html",
            s1="{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/read1_paired.fq_fastqc/summary.txt",
            s2="{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/read2_paired.fq_fastqc/summary.txt"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/fq.benchmark"
        shell:
            "fastqc {input.r1} {input.r2} --extract -o {params}"

    rule validateQCTrimm:
        input:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/read1_paired.fq_fastqc/summary.txt",
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/read2_paired.fq_fastqc/summary.txt",
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/read1_paired.fq_fastqc.html",
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/read1_paired.fq_fastqc.html",
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/fq_fw_internal_validation.txt",
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/fq_rv_internal_validation.txt"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/qc_validation.benchmark"
        script:
            "Scripts/validateQC.py"
else:
    rule skip_trimm_fastqc:
        output:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/fq_fw_internal_validation.txt",
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/fq_rv_internal_validation.txt"
        shell:
            "echo \"User skip raw reads quality control\" > {output[0]} && echo \"User skip raw reads quality control\" > {output[1]}"

if config["TAXONOMY"]["PROFILING"] == "KRAKEN" or config["TAXONOMY"]["PROFILING"] == "ALL":
    rule kraken:
        """
            Execute Kraken taxonomy profiling
        """
        input:
            fw="{PROJECT}/samples/{sample}/rawdata/fw.fastq"
            if config["TAXONOMY"]["KRAKEN"]["raw_reads"] == "Y" else "{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
            rv="{PROJECT}/samples/{sample}/rawdata/rv.fastq"
            if config["TAXONOMY"]["KRAKEN"]["raw_reads"] == "Y" else "{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.out"
        shell:
            "kraken --preload --db {config[TAXONOMY][KRAKEN][db]} --paired {input.fw} {input.rv} "
            "--threads {config[TAXONOMY][KRAKEN][threads]} {config[TAXONOMY][KRAKEN][extra_params]} > {output}"
    rule prepare_kraken_report:
        """
            Prepare the input file for addTaxonNames script
        """
        input:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.out"
        output:
            temp("{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.out.tmp")
        shell:
            "cat {input} | cut -f1,2,3 > {output}"
    rule kraken_labels:
        """
            Prepare the input file for addTaxonNames script
        """
        input:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.out.tmp"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.out.labels"
        shell:
            "addTaxonNames -t {config[TAXONOMY][KRAKEN][nodes]} -n {config[TAXONOMY][KRAKEN][names]} "
            "-i {input} {config[TAXONOMY][taxonomy_path]}  -o {output}"
    rule kraken_report:
        """
            Prepare the input file for addTaxonNames script
        """
        input:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.out.labels"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.report"
        shell:
            "kaijuReport -t {config[TAXONOMY][KRAKEN][nodes]} -n {config[TAXONOMY][KRAKEN][names]} "
            "-i {input} {config[TAXONOMY][taxonomy_path]}  -o {output}"

if config["TAXONOMY"]["PROFILING"] == "KAIJU" or config["TAXONOMY"]["PROFILING"] == "ALL":
    rule kaiju:
        """
            Execute Kaiju taxonomy profiling
        """
        input:
            fw="{PROJECT}/samples/{sample}/rawdata/fw.fastq"
            if config["TAXONOMY"]["KRAKEN"]["raw_reads"] == "Y" else "{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
            rv="{PROJECT}/samples/{sample}/rawdata/rv.fastq"
            if config["TAXONOMY"]["KRAKEN"]["raw_reads"] == "Y" else "{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kaiju.taxonomy.out"
        shell:
            "/export/data/aabdala/utils/kaiju/bin/kaiju -i {input.fw} -j {input.rv} "
            " -t {config[TAXONOMY][KAIJU][nodes]}  -f {config[TAXONOMY][KAIJU][db]} "
            "-z {config[TAXONOMY][KAIJU][threads]} {config[TAXONOMY][KAIJU][extra_params]} -o {output}"
    rule kaiju_labels:
        """
            addTaxonLabels for kaiju
        """
        input:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kaiju.taxonomy.out"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kaiju.taxonomy.out.labels"
        shell:
            "addTaxonNames -t {config[TAXONOMY][KAIJU][nodes]} -n {config[TAXONOMY][KAIJU][names]} "
            "-i {input} {config[TAXONOMY][taxonomy_path]}  -o {output}"
    rule kaiju_report:
        """
            addTaxonLabels for kaiju
        """
        input:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kaiju.taxonomy.out"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kaiju.taxonomy.out.labels"
        shell:
            "kaijuReport -t {config[TAXONOMY][KAIJU][nodes]} -n {config[TAXONOMY][KAIJU][names]} "
            "-i {input} {config[TAXONOMY][taxonomy_path]}  -o {output}"


# if config["TAXONOMY"]["PROFILING"] == "CLARK" or config["TAXONOMY"]["PROFILING"] == "ALL":
#     rule clark_target_db:
#         """
#             In order to run clark allways first we need to set the targets for the
#             database. For this we call the set_targets.sh script
#         """
#         output:
#             "{PROJECT}/runs/{run}/{sample}_data/taxonomy/clark.targetdb.out"
#         shell:
#             "set_targets.sh /export/data/databases/clark {config[TAXONOMY][CLARK][targets]} "
#             "{config[TAXONOMY][CLARK][level]} > {output}"
#     if config["TAXONOMY"]["CLARK"]["spaced"] == "--spaced":
#         rule clark_build_db:
#             """
#                 Clark can run in different modes, and one of this is the spaced db
#                 search. It requires more time and memory but it suppose it outperform
#                 the other taxonomy profielers, to run this spaced kmaer search it is
#                 needed to always run this command before run the classifier.
#             """
#             input:
#                 "{PROJECT}/runs/{run}/{sample}_data/taxonomy/clark.targetdb.out"
#             output:
#                 "{PROJECT}/runs/{run}/{sample}_data/taxonomy/clark.builddb.out"
#             shell:
#                 "buildSpacedDB.sh > {output}"
#     rule clark:
#         """
#             Execute Clark taxonomy profiling
#         """
#         input:
#             dbcreation_log="{PROJECT}/runs/{run}/{sample}_data/taxonomy/clark.builddb.out"
#             if config["TAXONOMY"]["CLARK"]["spaced"] == "--spaced" else "{PROJECT}/runs/{run}/{sample}_data/taxonomy/clark.targetdb.out",
#             fw="{PROJECT}/samples/{sample}/rawdata/fw.fastq"
#             if config["TAXONOMY"]["KRAKEN"]["raw_reads"] == "Y" else "{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
#             rv="{PROJECT}/samples/{sample}/rawdata/rv.fastq"
#             if config["TAXONOMY"]["KRAKEN"]["raw_reads"] == "Y" else "{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq"
#         params:
#             "{PROJECT}/runs/{run}/{sample}_data/taxonomy/"
#         output:
#             "{PROJECT}/runs/{run}/{sample}_data/taxonomy/clark.taxonomy.out"
#         shell:
#             "classify_metagenome.sh -P {input.fw} {input.rv} "
#             "{config[TAXONOMY][CLARK][extra_params]} -R {output} {config[TAXONOMY][CLARK][spaced]}"
if config["TAXONOMY"]["PROFILING"] not in "KRAKEN KAIJU CLARK ALL":
    rule create_taxo_out:
        """
            As there is no taxonomy profiling touch one file to generate a "silly" file
        """
        output:
            "{PROJECT}/runs/{run}/{sample}_data/no_tax.txt"
        shell:
            "touch {output}"
elif config["TAXONOMY"]["PROFILING"] == "ALL":
    rule merge_taxonomy_outs:
        """
            Merge the output of the three taxonomic profiler into one single outfile
        """
        input:
            kraken="{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.out",
            kaiju="{PROJECT}/runs/{run}/{sample}_data/taxonomy/kaiju.taxonomy.out"
            #clark="{PROJECT}/runs/{run}/{sample}_data/taxonomy/clark.taxonomy.out"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/taxonomy/all.taxonomy.out"
        shell:
            "touch {output}"

if config["ASSEMBLER"] == "SPADES":
    rule fq2fasta:
        input:
            read1_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
            read2_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq",
        output:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/reads_merged.fasta" if config["gzip_input"] == "F"
            else "{PROJECT}/runs/{run}/{sample}_data/trimmed/reads_merged.fastq.gz"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/fq2fasta.benchmark"
        shell:
            "fq2fa --merge {input.read1_paired} {input.read2_paired} {output}"

    rule concat_single_reads:
        input:
            read1_single="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_singles.fq",
            read2_single="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_singles.fq",
            tmp_flw="{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/fq_rv_internal_validation.txt"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/all_singles.fq"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/concat_single_reads.benchmark"
        shell:
            "cat {input.read1_single} {input.read2_single} > {output}"

    rule meta_spades:
        input:
            reads_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/reads_merged.fasta" if config["gzip_input"] == "F"
            else "{PROJECT}/runs/{run}/{sample}_data/trimmed/reads_merged.fastq.gz",
            read12_singles="{PROJECT}/runs/{run}/{sample}_data/trimmed/all_singles.fq"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta_tmp",
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta_tmp"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/assembly.benchmark"
        shell:
            "nice -{config[spades][nice]} spades.py --meta -t {config[spades][threads]} -m {config[spades][memory]} "
            "-k {config[spades][kmers]} --12 {input.reads_paired} -s {input.read12_singles} "
            "{config[spades][extra_params]} -o {params}"

    rule meta_spades_old:
        input:
            read1_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
            read2_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq",
            read12_singles="{PROJECT}/runs/{run}/{sample}_data/trimmed/all_singles.fq"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta",
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/assembly.benchmark"
        shell:
            "nice -{config[spades][nice]} spades.py --meta -t {config[spades][threads]} -m {config[spades][memory]} "
            "-k {config[spades][kmers]} --pe1-1 {input.read1_paired} --pe1-2 {input.read2_paired} --pe1-s {input.read12_singles} "
            "{config[spades][extra_params]} -o {params}"

if config["ASSEMBLER"] == "MEGAHIT":
    rule megahit:
        input:
            read1_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
            read2_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq",
            tmp_flw="{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/fq_rv_internal_validation.txt"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/final.contigs.fa"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]
        shell:
            "megahit -1 {input.read1_paired} -2 {input.read2_paired} -f --k-min {config[megahit][kmin]} "
            "--k-max {config[megahit][kmax]} --k-step {config[megahit][kstep]} {config[megahit][extra_params]} "
            "-m {config[megahit][memory]} -t {config[megahit][cpus]} -o {params}"
    rule std_assembly_megahit:
        input:
            contig="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/final.contigs.fa"
        output:
            contigs="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta",
            scaffolds="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
        shell:
            "mv {input.contig} {output.contigs} && ln -sr {output.contigs} {output.scaffolds}"
    #rule std_assembly_megahit_step2:
    #    input:
    #        contig="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta"
    #    output:
    #        scaffolds="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
    #    shell:
    #        "ln -sr {input.contig} scaffolds.fasta"

if config["ASSEMBLER"] == "IDBA":
    rule fq2fasta:
        input:
            read1_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
            read2_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq",
            tmp_flw="{PROJECT}/runs/{run}/{sample}_data/trimmed/qc/fq_rv_internal_validation.txt"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/reads_merged.fasta"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/fq2fasta.benchmark"
        shell:
            "fq2fa --merge {input.read1_paired} {input.read2_paired} {output}"
    #IN order to run idba it is needed to make somechanges into the source code:
    #https://groups.google.com/forum/#!topic/hku-idba/NE2JXqNvTFY and
    #http://seqanswers.com/forums/showthread.php?t=29109
    #change
    #static const uint32_t kMaxShortSequence = 128;
    #by
    #static const uint32_t kMaxShortSequence = 256;
    #For this reason, the pipe line uses my local installation, try to make it for all
    rule idba:
        input:
            "{PROJECT}/runs/{run}/{sample}_data/trimmed/reads_merged.fasta"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contig.fa",
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffold.fa"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/assembly.benchmark"
        shell:
            "/export/data/aabdala/utils/idba/bin/idba_ud -r {input} -o {params} "
            "--step {config[idba][step]} --num_threads {config[idba][threads]} {config[idba][extra_params]}"
    #As spades output contigs.fasta and scaffolds.fasta, we standarize those
    #names in order to decress complexity on downstream rules
    rule std_assembly:
        input:
            contig="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contig.fa",
            scaffold="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffold.fa"
        output:
            contigs="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta",
            scaffolds="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
        shell:
            "mv {input.contig} {output.contigs} && mv {input.scaffold} {output.scaffolds}"


rule quast_contigs:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/contigs/report.txt"
    params:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/contigs/"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/contigs/quast_contigs.benchmark"
    shell:
        "quast.py -t {config[quast][threads]} -o {params} {config[quast][extra_params]} {input}"

rule quast_scaffolds:
    input:
        scaffolds="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/scaffolds/report.txt"
    params:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/scaffolds/"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/scaffolds/quast_scaffolds.benchmark"
    shell:
        "quast.py -t {config[quast][threads]} -o {params} -s {input} --scaffolds {config[quast][extra_params]}"


if config["SPLIT_ASSEMBLY"] == "T":
    """
    This option will split the assembled reads into smaller chuncks.
    In order to do not affect the downstream rules (previously implemented)
    the splitted contig file is renamed as the original contig file and the original
    is renamed as contigs.complete.fasta at the end
    """
    rule split_assembly:
        input:
            contigs="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.chunks.fasta"
        shell:
            "/opt/biolinux/anaconda2.2019.07/bin/cut_up_fasta.py -c {config[SPLIT_SIZE]} "
            "-o 0 -m {input.contigs} > {output}"
    rule std_splitted_assembly:
        input:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.chunks.fasta"
        output:
            flag="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/split_flag.txt",
            complete="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.complete.fasta"
        params:
            contigs="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta"
        shell:
            "mv {params.contigs} {output.complete} && mv {input} {params.contigs} "
            "&& echo \"contigs.fasta has been splitted by cut_up_fasta.py original contig file is: contigs.complete.fasta\" > {output.flag}"
else:
    rule skip_split_assembly:
        output:
            flag="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/split_flag.txt"
        shell:
            "echo \"contigs.fasta has not been splitted\" > {output.flag}"

rule validate_assembly:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/scaffolds/report.txt",
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/contigs/report.txt",
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/split_flag.txt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/validate_assembly.txt"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/validate_assembly.benchmark"
    params:
        "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/"
    script:
        "Scripts/validateQuast.py"

if config["GENE_CALLING"]["TOOL"] == "MGM":
    #rule check_mgm_key:
    #    output:
    #        "$HOME/.gm_key"
    #    shell:
    #        "cp {config[GENE_CALLING][MGM][key]}  {output}"
    rule call_genes_mgm:
        input:
            assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta",
            #key="$HOME/.gm_key"
        output:
            gff="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/mgm_"+config["ANALYSIS"]+"/genes.gff",
            genes="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/mgm_"+config["ANALYSIS"]+"/genes.fasta",
            prots="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/mgm_"+config["ANALYSIS"]+"/prots.fasta"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/mgm_"+config["ANALYSIS"]+"/"
        shell:
            "/export/data/aabdala/utils/MetaGeneMark_linux_64/mgm/gmhmmp -a -d -f {config[GENE_CALLING][MGM][f]} "
            "-m {config[GENE_CALLING][MGM][model]} -A {output.genes} -D {output.prots} -o {output.gff} {input.assembly}"

elif config["GENE_CALLING"]["TOOL"] == "FGS":
    rule call_genes_fgs:
        input:
            assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta"
        output:
            gff="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/fgs_"+config["ANALYSIS"]+"/genes.gff",
            genes="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/fgs_"+config["ANALYSIS"]+"/genes.ffn",
            prots="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/fgs_"+config["ANALYSIS"]+"/genes.faa",
            out="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/fgs_"+config["ANALYSIS"]+"/genes.out"


        params:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/fgs_"+config["ANALYSIS"]+"/"
        shell:
            "/export/data/aabdala/utils/FragGeneScan/run_FragGeneScan.pl -genome={input.assembly} -out={params}genes "
            "-complete={config[GENE_CALLING][FGS][complete]} -train={config[GENE_CALLING][FGS][model]} -thread={config[GENE_CALLING][FGS][threads]}"

elif config["GENE_CALLING"]["TOOL"] == "PRODIGAL":
    rule call_genes_prodigal:
        input:
            assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta"
        output:
            gff="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/prodigal_"+config["ANALYSIS"]+"/genes.gff",
            genes="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/prodigal_"+config["ANALYSIS"]+"/genes.ffn",
            prots="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/prodigal_"+config["ANALYSIS"]+"/genes.faa"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/prodigal_"+config["ANALYSIS"]+"/"
        shell:
            "prodigal -a {output.prots} -d {output.genes} -f gff  -i {input.assembly} -p meta -o {output.gff} "
else :
    rule skip_gene_calling:
        output:
             "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/gene_calling.skip"
        shell:
            "touch {output}"

  #cp gm_key_64 ~/.gm_key
rule bwa_index:
    input:
        tmp_flw="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/quast/validate_assembly.txt",
        assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
        if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_assembly.bwt"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/bwaindex.benchmark"
    params:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_assembly"
    shell:
        "bwa index -p {params} {input.assembly}"
#nice -5 bwa mem -t 4 $folderOut/cross-assembly $fq1 $fq2 \
#   | samtools view -b - \
#   | samtools sort - -o $fol/mapped_against_cross-assembly_sorted.bam
rule bwa_mem:
    input:
        read1_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read1_paired.fq",
        read2_paired="{PROJECT}/runs/{run}/{sample}_data/trimmed/read2_paired.fq",
        tmp_flw="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_assembly.bwt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_mapped_against_cross-assembly_sorted.bam"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/bwamem.benchmark"
    params:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_assembly"
    shell:
        "nice -{config[bwa][nice]} bwa mem -t {config[bwa][threads]} {params} "
        "{input.read1_paired} {input.read2_paired} "
        "| samtools view -b - | samtools sort - -o {output}"

rule sam_flags:
     input:
         "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_mapped_against_cross-assembly_sorted.bam"
     output:
         "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_mapped_against_cross-assembly_sorted.flagstat"
     benchmark:
         "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/samflags.benchmark"
     shell:
         "samtools flagstat {input} > {output}"
rule summarize_bam:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_mapped_against_cross-assembly_sorted.bam"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/jgi_summ.benchmark"
    shell:
        "jgi_summarize_bam_contig_depths --outputDepth {output} {config[jgi_summ][extra_params]} {input}"
rule avg_coverage:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth_avg.txt"
    shell:
        "cut -f1,3 {input} > {output}"

if config["BINNING"] == "METABAT" or config["BINNING"] == "DAS":
    rule metabat:
        input:
            depth="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt",
            assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta"
        output:
        #"{PROJECT}/runs/{run}/{sample}_data/metabat2/bin/metabat2.log",
        #"{PROJECT}/runs/{run}/{sample}_data/metabat2/bin{number}.fa",
        #expand("{PROJECT}/runs/{run}/{sample}_data/metabat2/bin{{number}}.fa",PROJECT=config["PROJECT"],sample=config["SAMPLES"], run=run),
        #dynamic("{PROJECT}/runs/{run}/{sample}_data/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/bin.{number}.fa")
        #"{PROJECT}/runs/{run}/{sample}_data/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/bin.{number}.fa"
            "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/metabat.log"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/bin"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/metabat.benchmark"
        shell:
            "metabat2 -o {params} -i {input.assembly} -t {config[metabat2][threads]} "
            "-m {config[metabat2][min_contig]} -a {input.depth} "
            "--minS {config[metabat2][min_score]} --maxP {config[metabat2][maxP]} "
            "--maxEdges {config[metabat2][maxEdge]} -s {config[metabat2][min_bin_size]} "
            "{config[metabat2][extra_params]}  > {output}"
if config["BINNING"] == "MAXBIN" or config["BINNING"] == "DAS":
    rule maxbin:
        """
        Please make sure that your abundance information is provided in the following format:
        (contig header)\t(abundance)
        """
        input:
            depth="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth_avg.txt",
            assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta"
        output:
            log="{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/maxbin.log"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/bin"
        benchmark:
            "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/maxbin.benchmark"
        shell:
            "/export/data/aabdala/utils/MaxBin-2.2.4/run_MaxBin.pl -contig {input.assembly} "
            #"run_MaxBin.pl -contig {input.assembly} "
            "-abund  {input.depth} -out {params} -thread {config[maxbin][threads]} "
            "-prob_threshold {config[maxbin][prob_threshold]} -markerset {config[maxbin][markerset]} "
            "-min_contig_length {config[maxbin][min_contig_length]}  "
            "{config[maxbin][plotmarker]} {config[maxbin][extra_params]} > {output.log}"
if config["BINNING"] == "CONCOCT" or config["BINNING"] == "DAS":
    """
    This rules try to follow the steps recommended by using CONCOCT according to
    the following documentation: https://concoct.readthedocs.io/en/latest/complete_example.html
    """
    rule concoct:
        input:
            depth="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth_avg.txt",
            assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta"
        output:
            #log="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/concoct.log",
            clustering="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/bin_clustering_gt"+config["concoct"]["min_contig_length"]+".csv"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/bin"
        shell:
            "concoct -l {config[concoct][min_contig_length]} -i {config[concoct][max_iteration]} "
            "--coverage_file {input.depth} --composition_file {input.assembly} -b {params} "
    """
    We only need the bins in case that the user is running concoct alone, otherwise
    we only use the clustering file to DAS to create new bins
    """
    if config["BINNING"] == "CONCOCT":
        rule extract_concoct_bins:
            input:
                assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
                if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta",
                clustering="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/bin_clustering_gt"+config["concoct"]["min_contig_length"]+".csv"
            output:
                log="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/concoct.log"
            params:
                "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/"
            shell:
                #"/usr/local/bioinf/concoct-0.4/scripts/extract_fasta_bins.py  --output_path  {params} {input.assembly} {input.clustering} > {output.log}"
                #"/export/data/aabdala/utils/CONCOCT/scripts/extract_fasta_bins.py  --output_path  {params} {input.assembly} {input.clustering} > {output.log}"
                "python /opt/biolinux/anaconda2.2019.07/pkgs/concoct-1.1.0-py27h88e4a8a_0/bin/extract_fasta_bins.py  --output_path  {params} {input.assembly} {input.clustering} > {output.log}"
if config["BINNING"] == "BINSANITY" or config["BINNING"] == "DAS":
    rule log_transform_coverage:
        input:
            depth="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth.txt"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth_avg_log.txt"
        shell:
            "cat {input} | awk 'NR>1 && $3>1{{ if($3 <= 1) a = 0; else  a = log($3)/log(10); printf(\"%s\\t%0.4f\\n\",$1,a)}}' > {output}"
            #"cat {input} | awk 'NR>1 {{print $1,log($3+1)/log(10)}}' > {output}"

    rule filter_fasta_by_coverage:
        input:
            depth="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth_avg_log.txt",
            fasta="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta"
        output:
            temp("{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/assembly_coverage_gt0.fasta")
        shell:
            "cat {input.depth} | cut -f1 | grep -w -A1 --no-group-separator -F -f - {input.fasta}  > {output}" 

    rule bin_sanity:
        input:
            depth="{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_depth_avg_log.txt",
            assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta"
        params:
            contig_directory="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"],
            bin_directory="{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/",
            #final bins will be at  {bin_directory}/BinSanity-Final-bins/final_Bin-xx.fna
            #also they can be found in the form final_Bin-xx_refined-xx.fna
            assembly="scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "contigs.fasta"
        output:
            log="{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanityWf.log"
        shell:
            "/export/data/aabdala/utils/BInSanity/BinSanity-master/bin/Binsanity-wf -f {params.contig_directory} -l {params.assembly} -c {input.depth} -o {params.bin_directory} --binPrefix  final "
            "-p {config[binsanity][preference]} -x {config[binsanity][min_contig_length]} --threads {config[binsanity][threads]} "
            "{config[binsanity][extra_params]}"
if config["BINNING"] == "DAS":
    rule bins_to_table:
        input:
            maxbin="{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/maxbin.log",
            metabat="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/metabat.log",
            #concoct="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/concoct.log"
            concoct="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/bin_clustering_gt"+config["concoct"]["min_contig_length"]+".csv",
            binsanity="{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanityWf.log"
        output:
            metabat_out="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv",
            maxbin_out="{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv",
            concoct_out="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv",
            binsanity_out="{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv"
        params:
            output_dir_metabat="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"],
            output_dir_maxbin="{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"],
            output_dir_concoct="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"],
            output_dir_binsanity="{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanity-Final-bins/",
            file_ext_metabat="fa",
            file_ext_maxbin="fasta",
            file_ext_concoct="fa",
            file_ext_binsanity="fna"
        script:
            "Scripts/tableBins.py"
    rule das:
        input:
            metabat_bin2t="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv",
            maxbin_bin2t="{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv",
            concoct_bin2t="{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv",
            binsanity_bin2t="{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/binTable.tsv",
            assembly="{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/scaffolds.fasta"
            if config["ANALYSIS"] == "SCAFFOLDS" else "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/contigs.fasta"
        params:
            "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/DasOut"
        output:
            "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log"
        shell:
            "DAS_Tool -i {input.metabat_bin2t},{input.maxbin_bin2t},{input.concoct_bin2t},{input.binsanity_bin2t} "
            "-l metabat,maxbin,concoct,binsanity -c {input.assembly} -t {config[das][threads]} "
            "--write_bins 1 --db_directory {config[das][db]} --search_engine {config[das][search_engine]} "
            "--create_plots 1 {config[das][extra_params]} -o {params} > {output}"
else:
    rule skip_das:
        params:
            "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]
        output:
            "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log"
        shell:
            "touch {output}"

rule checkM_bins:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/metabat.log"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/maxbin.log"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/concoct.log"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanityWf.log"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log"
    params  :
        bin_folder="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanity-Final-bins/"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/DasOut_DASTool_bins/",
        bin_ext="fa"
        if config["BINNING"] == "METABAT" or config["BINNING"] == "CONCOCT" or config["BINNING"] == "DAS" else
        "fasta" if config["BINNING"] == "MAXBIN" else "fna",
        out_folder="{PROJECT}/runs/{run}/{sample}_data/binning/checkM/"
    output:
        out_file="{PROJECT}/runs/{run}/{sample}_data/binning/checkM/summary.txt"
    shell:
        "checkm lineage_wf -f {output.out_file} -t  {config[checkM][threads]} -x {params.bin_ext} {config[checkM][extra_params]} {params.bin_folder} {params.out_folder} "


rule prokka_bins:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/metabat.log"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/maxbin.log"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/concoct.log"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanityWf.log"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log",
        "{PROJECT}/runs/{run}/{sample}_data/binning/checkM/summary.txt"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
    params:
        output_dir="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanity-Final-bins/"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/DasOut_DASTool_bins/",
        file_ext= "fa"
        if config["BINNING"] == "METABAT" or config["BINNING"] == "CONCOCT" or config["BINNING"] == "DAS" else
        "fasta"
        if config["BINNING"] == "BINSANITY" else
        "fna"
    benchmark:
        #"{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.benchmark"
        #if config["BINNING"] == "METABAT" else
        #"{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.benchmark"
        "{PROJECT}/runs/{run}/{sample}_data/binning/prokka"+config["BINNING"]+"_bins.benchmark"
    script:
            "Scripts/annotateProkkaBins.py"
rule diamond_bins:
    input:
        #"{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        #if config["BINNING"] == "METABAT" else
        #"{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
    params:
        output_dir="{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/BinSanity-Final-bins/"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/DasOut_DASTool_bins/",
        file_ext= "faa"
    benchmark:
        "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.benchmark"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.benchmark"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.benchmark"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.benchmark"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.benchmark"
    script:
        "Scripts/diamondProkkaBins.py"

rule report:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"_mapped_against_cross-assembly_sorted.flagstat",
        #"{PROJECT}/runs/{run}/{sample}_data/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/prokka.log",
        "{PROJECT}/runs/{run}/{sample}_data/binning/metabat2/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
        if config["BINNING"] == "METABAT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
        if config["BINNING"] == "MAXBIN" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/concoct/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
        if config["BINNING"] == "CONCOCT" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/binsanity/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log"
        if config["BINNING"] == "BINSANITY" else
        "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/diamond.log",
        #"{PROJECT}/runs/{run}/{sample}_data/binning/maxbin/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/maxbin.log",
        "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kraken.taxonomy.report"
         if config["TAXONOMY"]["PROFILING"] == "KRAKEN" else
         "{PROJECT}/runs/{run}/{sample}_data/taxonomy/kaiju.taxonomy.out.labels"
         if config["TAXONOMY"]["PROFILING"] == "KAIJU" else
         "{PROJECT}/runs/{run}/{sample}_data/taxonomy/clark.taxonomy.out"
         if config["TAXONOMY"]["PROFILING"] == "CLARK" else
         "{PROJECT}/runs/{run}/{sample}_data/taxonomy/all.taxonomy.out"
         if config["TAXONOMY"]["PROFILING"] == "ALL" else
         "{PROJECT}/runs/{run}/{sample}_data/no_tax.txt",
         "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/mgm_"+config["ANALYSIS"]+"/prots.fasta"
         if config["GENE_CALLING"]["TOOL"] == "MGM" else
         "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/fgs_"+config["ANALYSIS"]+"/genes.out"
         if config["GENE_CALLING"]["TOOL"] == "FGS" else
         "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/prodigal_"+config["ANALYSIS"]+"/genes.gff"
         if config["GENE_CALLING"]["TOOL"] == "PRODIGAL" else
         "{PROJECT}/runs/{run}/{sample}_data/assembly_"+config["ASSEMBLER"]+"/gene_calling.skip",
         "{PROJECT}/runs/{run}/{sample}_data/binning/das/"+config["ANALYSIS"]+"_"+config["ASSEMBLER"]+"/das.log"
        #"{PROJECT}/runs/{run}/{sample}_data/metabat2/prokka.out"
        #"{PROJECT}/runs/{run}/{sample}_data/metabat2/bin/metabat2.log"
    output:
        temp("{PROJECT}/runs/{run}/{sample}_data/report.html")
    script:
        "Scripts/report.py"
rule tune_report:
    input:
        "{PROJECT}/runs/{run}/{sample}_data/report.html"
    output:
        "{PROJECT}/runs/{run}/{sample}_data/report_f.html"
    script:
        "Scripts/tuneReport.py"
