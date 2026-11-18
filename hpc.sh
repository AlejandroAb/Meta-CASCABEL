#!/bin/sh
#SBATCH --partition=HMEM
#SBATCH --nodes=1
#SBATCH --mail-user=                                # Fill in your own email
#SBATCH --mail-type=FAIL,BEGIN,END,TIME_LIMIT_80
#SBATCH --time=5-10:00:00
#SBATCH --cpus-per-task=100                         # Keep this on 100 to use a full node and avoid crashing
#SBATCH --mem=1000G                                 # Maximum for HMEM

module load anaconda/2024.02
conda activate snakemake_v7.14.2
export GTDBTK_DATA_PATH="/export/lv13/databases/gtdb/release232"

# Run command:
snakemake --configfile config.yaml -c 100 --use-conda --conda-frontend conda 
# For a dry run, add to the above command: -np
# To continue a stopped run without starting all over, add to the above command: --rerun-triggers mtime

# Generate report file, 'report.hpc.zip' can be changed to any name you want
snakemake --configfile config.yaml --report report.hpc.zip

# These lines remove all bam and depth files at the end of the run:
rm -rf "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/*.bam"
rm -rf "{PROJECT}/runs/{run}/{sample}_data/bwa-mem/*.txt"
# If you want to keep these, comment them out by adding '#'' at the beginning of the line
# Also know that at the end of the run a file 'contig_coverage.txt' is created per sample in the directory: '{SAMPLE}_data/binning/FinalBins'
# Which contains the avarage depth of all contigs within the final bins
