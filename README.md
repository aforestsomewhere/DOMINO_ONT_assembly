# Snakemake workflow: DOMINO_ONT_assembly

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/aforestsomewhere/DOMINO_ONT_assembly/workflows/Tests/badge.svg?branch=main)](https://github.com/aforestsomewhere/DOMINO_ONT_assembly/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for the assembly of bacterial genomes sequenced on ONT platforms

## Usage

1. Prepare a working directory, containing a folder "fastq" which contains one fastq.gz file for each barcode, the "sequencing_summary_XXXXXX".txt file from the run (the config file will detect the full name) and a genome_sizes.csv file. This file should be comma separated, where the first column contains the barcodes, second column contains any isolate references, and the third column contains the expected genome size (in the format X.Xy, e.g. 4.9m for a 4.9Mb genome)
2. Download and make executable the "generate_config.sh" file
3. Run (./generate_config.sh) to generate the config file
4. Ensure the required conda environments have been set up.
5. source activate snakemake_env
6. snakemake --snakefile Snakefile.py --configfile config.yaml --jobs 40 --use-conda --debug-dag --cluster "qsub -V -cwd -pe thread {threads} -N {rule} -q {params.queue}" all
7. conda deactivate

# TODO

1. Comment!
2. Upgrade snakemake to version that supports directory() and checkpoints
