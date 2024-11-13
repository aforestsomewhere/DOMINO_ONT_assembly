# Snakemake workflow: DOMINO_ONT_assembly

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/aforestsomewhere/DOMINO_ONT_assembly/workflows/Tests/badge.svg?branch=main)](https://github.com/aforestsomewhere/DOMINO_ONT_assembly/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for the assembly of bacterial genomes sequenced on ONT platforms

## Usage

1. Prepare your reads by gathering ONT data in the ont_fastq folder, and Illumina data in the raw_fastq folder
2. Prepare your genome_sizes.csv file. A sample file is provided. This file should be comma separated, where the first column contains the barcodes, second column contains any isolate references, and the third column contains the expected genome size (in the format X.Xy, e.g. 4.9m for a 4.9Mb genome)
4. Run (./generate_config.sh) to generate the config file. Provide the flag --ont for ont data only, --illumina for illumina data only or --ont --illumina for both. Other variables, such as shell executable, conda profile and database paths may be modified within this script or directly within the config file (config/config.yaml)
5. Ensure the required conda environments have been set up for each tool ("porechop_abi_env", "filtlong_env", "nanocomp_env", "flye_2.9")
6. source activate snakemake_env
7. To run the QC pipeline - adaptor removal with Porechop, light QC with Filtlong, readset analysis with NanoPlot and NanoComp, and an initial first pass assembly with Flye: snakemake --snakefile workflow/ONT_QC.smk --config.yaml --profile smk-profile-slurm/ --use-conda --executor cluster-generic
8. conda deactivate

# TODO
1. Upgrade snakemake to version that supports directory() and checkpoints
