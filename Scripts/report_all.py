import subprocess
from snakemake.utils import report

def readBenchmark( benchFile ):
   with open(benchFile) as bfile:
       txt = "**Benchmark info:**\n\n"
       for l in bfile:
         txt += l + "\n\n"
   return txt;

################
#Function to retrive the sample names and put in the report title
#@param file with the sample list, it is created during combine_filtered_samples
#snakemake.wildcards.project + "/runs/" + snakemake.wildcards.run + "/samples.log"
#@return the title with the samples
def getSampleList(sampleFile):
    with open(sampleFile) as sfile:
        samps ="Amplicon Analysis Report for Libraries: "
        for l in sfile:
            samps+= l
        samps+="\n"
        for i in range(0,len(samps)):
            samps+="="
    return samps;

def getCombinedSamplesList(sampleFile):
    with open(sampleFile) as sfile:
        command =":commd:`"
        i=0
        for l in sfile:
            if i == 0:
                command+= l + "`\n\n"
            i+=1
    return command;

#fqv = subprocess.run(['fastqc', '--version'], stdout=subprocess.PIPE)
#fqVersion = "**" + fqv.stdout.decode('utf-8').strip() + "**"

title = getSampleList(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/samples.log")
catCommand =  getCombinedSamplesList(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/cat_samples.log")

totalReads = "TBD"
#ps = subprocess.Popen(("grep", "-e '^>' "+snakemake.wildcards.PROJECT+"/samples/"+snakemake.wildcards.run+"/seqs_fw_rev_filtered.fasta"), stdout=subprocess.PIPE)
#totalReads = subprocess.check_output(('wc', '-l'), stdin=ps.stdout)
#ps.wait()

combineBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/combine_seqs_fw_rev.benchmark")
otuBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu.benchmark")
pikRepBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/pick_reps.benchmark")
assignTaxaBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/assign_taxa.benchmark")
otuTableBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/otuTable.biom.benchmark")
convertOtuBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/otuTable.txt.benchmark")
otuNoSingletonsBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/otuTable_nosingletons.bio.benchmark")
#filter_fasta.py -f {input.fastaRep} -o {output} -b {input.otuNoSingleton} {config[filterFasta][extra_params]}
filterBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/representative_seq_set_noSingletons.benchmark")
#align_seqs.py -m {config[alignRep][m]} -i {input} -o {params.outdir} {config[alignRep][extra_params]}
alignBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/aligned/align_rep_seqs.benchmark")
#"filter_alignment.py -i {input} -o {params.outdir} {config[filterAlignment][extra_params]}"
alignFilteredBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/aligned/filtered/align_rep_seqs.benchmark")
#"make_phylogeny.py -i {input} -o {output} {config[makeTree][extra_params]}"
makePhyloBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/aligned/filtered/representative_seq_set_noSingletons_aligned_pfiltered.benchmark")
#get versions
clusterOtuV = subprocess.run(['pick_otus.py', '--version'], stdout=subprocess.PIPE)
clusterOtuVersion = "**" + clusterOtuV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

pickRepV = subprocess.run(['pick_rep_set.py', '--version'], stdout=subprocess.PIPE)
pickRepVersion = "**" + pickRepV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

assignTaxaV = subprocess.run(['parallel_assign_taxonomy_blast.py', '--version'], stdout=subprocess.PIPE)
assignTaxaVersion = "**" + assignTaxaV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

makeOTUV = subprocess.run(['make_otu_table.py', '--version'], stdout=subprocess.PIPE)
makeOTUVersion = "**" + makeOTUV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

convertBiomV = subprocess.run(['biom', '--version'], stdout=subprocess.PIPE)
convertBiomVersion = "**" + convertBiomV.stdout.decode('utf-8').strip() + "**"

filterOTUNoSV = subprocess.run(['filter_otus_from_otu_table.py', '--version'], stdout=subprocess.PIPE)
filterOTUNoSVersion = "**" + filterOTUNoSV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

filterFastaV = subprocess.run(['filter_fasta.py', '--version'], stdout=subprocess.PIPE)
filterFastaVersion = "**" + filterFastaV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

alignFastaVersion="TBD"
try:
    alignFastaV = subprocess.run(['align_seqs.py', '--version'], stdout=subprocess.PIPE)
    if "Version" in alignFastaVersion:
        alignFastaVersion = "**" + alignFastaV.stdout.decode('utf-8').replace('Version: ','').strip() + "**"
except Exception as e:
    alignFastaVersion = "**Problem retriving the version**"

filterAlignmentV = subprocess.run(['filter_alignment.py', '--version'], stdout=subprocess.PIPE)
filterAlignmentVersion = "**" + filterAlignmentV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

makePhyloV = subprocess.run(['make_phylogeny.py', '--version'], stdout=subprocess.PIPE)
makePhyloVersion = "**" + makePhyloV.stdout.decode('utf-8').replace('Version:','').strip() + "**"
try:
    treads = subprocess.run( ["grep '^>' " + snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/seqs_fw_rev_filtered.fasta | wc -l"], stdout=subprocess.PIPE, shell=True)
    totalReads = "**" + treads.stdout.decode('utf-8').strip() + "**"
except Exception as e:
    totalReads = "Problem reading outputfile"

intOtus = 1
try:
    totus = subprocess.run( ["cat " +  snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/otu/seqs_fw_rev_filtered_otus.txt | wc -l"], stdout=subprocess.PIPE, shell=True)
    intOtus = int(totus.stdout.decode('utf-8').strip())
    #print("Total OTUS" + str(intOtus))
    totalOtus = "**" + str(intOtus) + "**"
except Exception as e:
    totalOtus = "**Problem reading outputfile**"

prcAssigned = 0.0
prcNotAssignedOtus="TBD"
try:
    aOtus = subprocess.run( ["grep 'No blast hit' " +  snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/otu/representative_seq_set_tax_assignments.txt | wc -l"], stdout=subprocess.PIPE, shell=True)
    notAssignedOtus = int(aOtus.stdout.decode('utf-8').strip())
    #print("Not assigned OTUS" + str(notAssignedOtus))
    prcAssigned = ((intOtus - notAssignedOtus)/intOtus)*100
    prcAssignedOtus = "**" + "{:.2f}".format(prcAssigned) + "%**"
except Exception as e:
    prcAssignedOtus = "**Problem reading outputfile**"


#include user description on the report
desc = snakemake.config["description"]
txtDescription = ""
if len(desc) > 0:
    txtDescription = "\n**User description:** "+desc+"\n"

#report benchmark
#reportBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/samples/"+snakemake.wildcards.run+"/report_all.benchmark")
#:commd:`pick_otus.py -m {snakemake.config[pickOTU][m]} -i {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/seqs_fw_rev_filtered.fasta -o {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/  -s {snakemake.config[pickOTU][s]} {snakemake.config[pickOTU][extra_params]}`

report("""
{title}
    .. role:: commd
    .. role:: red
    .. role:: green

This report consists on the OTU creation and taxonomic assignation for all the combined accepted reads of given samples.
{txtDescription}
The total amount of reads is: {totalReads}

Combine Reads
---------------

This is the first step to perform the taxonomic analysis, which consists on merge all the libraries reads into one single file.

This step was performed with the following command:

{catCommand}

**The output file for this step is:**

:green:`- Merged reads:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/seqs_fw_rev_filtered.fasta

{combineBenchmark}

Cluster OTUs
-------------

The OTU picking step assigns similar sequences to operational taxonomic units, or OTUs, by clustering sequences based on a user-defined similarity threshold. Sequences which are similar at or above the threshold level are taken to represent the presence of a taxonomic unit (e.g., a genus, when the similarity threshold is set at 0.94) in the sequence collection.

Currently, QIIME permits to perform this step with different clustering methods in this case, the one used was: {snakemake.config[pickOTU][m]}

:commd:`pick_otus.py -m {snakemake.config[pickOTU][m]} -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/seqs_fw_rev_filtered.fasta -o {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/  -s {snakemake.config[pickOTU][s]} {snakemake.config[pickOTU][extra_params]} .`

:red:`The version used for this run was:` {clusterOtuVersion}

**The output files for this step are:**

:green:`- OTU List:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/seqs_fw_rev_filtered_otus.txt

:green:`- Log file:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/seqs_fw_rev_filtered_otus.log

The total number of different OTUS is: {totalOtus}

{otuBenchmark}

Pick representative OTUs
--------------------------
After picking OTUs, it is necessary to pick a representative set of sequences. For each OTU, you will end up with one sequence that can be used in subsequent analyses.

:commd:`pick_rep_set.py -m {snakemake.config[pickRep][m]} -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/seqs_fw_rev_filtered_otus.txt -f {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/seqs_fw_rev_filtered.fasta -o {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/representative_seq_set.fasta --log_fp {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/representative_seq_set.log {snakemake.config[pickRep][extra_params]} .`

:red:`The version used for this run was:` {pickRepVersion}

**The output files for this step are:**

:green:`- Fasta file with representative sequences:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/representative_seq_set.fasta

:green:`- Log file:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/representative_seq_set.log

{pikRepBenchmark}

Assign taxonomy
----------------
Given a set of sequences, this step attempts to assign the taxonomy of each sequence. This script performs like the assign_taxonomy.py script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel.

:commd:`parallel_assign_taxonomy_blast.py -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/representative_seq_set.fasta --id_to_taxonomy_fp {snakemake.config[assignTax][mappFile]} --reference_seqs_fp {snakemake.config[assignTax][dbFile]} --jobs_to_start {snakemake.config[assignTax][jobs]} --output_dir {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/  {snakemake.config[assignTax][extra_params]}.`

:red:`The version used for this run was:` {assignTaxaVersion}

The percentage of successfully assigned OTUs is: {prcAssignedOtus}

**The output files for this step are:**

:green:`- OTU taxonomy assignation:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/representative_seq_set_tax_assignments.txt

:green:`- Log file:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/representative_seq_set_tax_assignments.log

{assignTaxaBenchmark}

Make OTU table
---------------
The script make_otu_table.py tabulates the number of times an OTU is found in each sample, and adds the taxonomic predictions for each OTU in the last column if a taxonomy file is supplied.

:commd:`make_otu_table.py -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/seqs_fw_rev_filtered_otus.txt -t {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/representative_seq_set_tax_assignments.txt -o {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/otuTable.biom {snakemake.config[makeOtu][extra_params]}.`

:red:`The version used for this run was:` {makeOTUVersion}

**The output file for this step is:**

:green:`- Biom format table:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/otuTable.biom

{otuTableBenchmark}

Convert OTU table
------------------
Convert from the BIOM table format to specific target format

:commd:`biom convert -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/otuTable.biom -o {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/otuTable.txt {snakemake.config[biom][tableType]} {snakemake.config[biom][headerKey]} {snakemake.config[biom][outFormat]} {snakemake.config[biom][extra_params]}.`

:red:`The version used for this run was:` {convertBiomVersion}

**The output file for this step is:**

:green:`- TSV format table:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/otuTable.txt

{convertOtuBenchmark}

Filter OTU table
-----------------
Filter OTUs from an OTU table based on their observation counts or identifier. So, in this case we remove all the OTUs that have less than {snakemake.config[filterOtu][n]} observation counts.

:commd:`filter_otus_from_otu_table.py -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/otuTable.biom -o {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/otuTable_noSingletons.biom -n {snakemake.config[filterOtu][n]} {snakemake.config[filterOtu][extra_params]}.`

:red:`The version used for this run was:` {filterOTUNoSVersion}

**The output file for this step is:**

:green:`- Biom format table:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/otuTable_noSingletons.biom*

{otuNoSingletonsBenchmark}

Convert Filtered OTU table
---------------------------
Convert the filtered OTU table from the BIOM table format to specific target format

:commd:`biom convert -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/otuTable_noSingletons.biom -o {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/otuTable_noSingletons.txt {snakemake.config[biom][tableType]} {snakemake.config[biom][headerKey]} {snakemake.config[biom][outFormat]} {snakemake.config[biom][extra_params]}.`

:red:`The version used for this run was:` {convertBiomVersion}

**The output file for this step is:**

:green:`- TSV format table:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/otuTable_noSingletons.txt

{otuNoSingletonsBenchmark}

Filter representative sequences
---------------------------------
This script can be applied to remove sequences from a fasta or fastq file based on input criteria. In this case, we retain the sequence according to the filtered OTU biom table

:commd:`filter_fasta.py -f {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/representative_seq_set.fasta -o {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/representative_seq_set_noSingletons.fasta -b {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/otuTable_noSingletons.biom {snakemake.config[filterFasta][extra_params]}.`

:red:`The version used for this run was:` {filterFastaVersion}

**The output file for this step is:**

:green:`- Filtered fasta file:` {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/representative_seq_set_noSingletons.fasta

{filterBenchmark}

Align representative sequences
-------------------------------
This script aligns the sequences in a FASTA file to each other or to a template sequence alignment, depending on the method chosen. In this case the selected method was: {snakemake.config[alignRep][m]}

:commd:`align_seqs.py -m {snakemake.config[alignRep][m]} -i  {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/representative_seq_set_noSingletons.fasta -o  {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/aligned/representative_seq_set_noSingletons_aligned.fasta {snakemake.config[alignRep][extra_params]}.`

:red:`The version used for this run was:` {alignFastaVersion}

**The output files for this step are:**

:green:`- Aligned fasta file:` {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/aligned/representative_seq_set_noSingletons_aligned.fasta
:green:`- Log file:` {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/aligned/representative_seq_set_noSingletons_log.txt

{alignBenchmark}

Filter alignment
-----------------
Filter sequence alignment by removing highly variable regions.

This script should be applied to generate a useful tree when aligning against a template alignment

:commd:`filter_alignment.py -i  {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/aligned/representative_seq_set_noSingletons_aligned.fasta -o {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/aligned/filtered/ {snakemake.config[filterAlignment][extra_params]} .`

:red:`The version used for this run was:` {filterAlignmentVersion}

**The output file for this step is:**

:green:`- Aligned fasta file:` {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/aligned/filtered/representative_seq_set_noSingletons_aligned_pfiltered.fasta

{alignFilteredBenchmark}

Make tree
-----------
Many downstream analyses require that the phylogenetic tree relating the OTUs in a study be present. The script make_phylogeny.py produces this tree from a multiple sequence alignment. Trees are constructed with a set of sequences representative of the OTUs, in this case with **{snakemake.config[makeTree][method]}** method.

:commd:`make_phylogeny.py -i {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/aligned/representative_seq_set_noSingletons_aligned_pfiltered.fasta -o {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/aligned/representative_seq_set_noSingletons_aligned_pfiltered.tre -t {snakemake.config[makeTree][method]} {snakemake.config[makeTree][extra_params]}.`

:red:`The version used for this run was:` {makePhyloVersion}

**The output file for this step is:**

:green:`- Taxonomy tree:` {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/aligned/representative_seq_set_noSingletons_aligned_pfiltered.tre

{makePhyloBenchmark}

""", snakemake.output[0], metadata="Author: J. Engelmann & A. Abdala ")
