# DOMINO WP3-4-Assemblies-ONT
# Assembly for low-depth isolates not fit for sub-sampling
# Author: Katie O'Mahony
# Adapted from: Julien Tap

shell.executable(config['shell_executable'])
shell.prefix(config['shell_prefix'])

import os
import glob
import warnings
import gzip

# fastq_path = config['raw_fastq_path']
genome_sizes = config['genome_sizes']
###HELPER FUNCTIONS #####

# Define the workflow
rule all:
    input:
        expand(os.path.join(config['assemblies_work_path_shallow'], "{sample}", "flye","assembly.fasta"),sample=config['shallow_sample_names']),
        expand(os.path.join(config['assemblies_work_path_shallow'], "{sample}", "raven","assembly.fasta"),sample=config['shallow_sample_names']),
        expand(os.path.join(config['assemblies_work_path_shallow'], "{sample}", "unicycler","assembly.fasta"),sample=config['shallow_sample_names'])

# Rule for assembling with flye using subsamples
rule flye_shallow:
    input:
        fastq=os.path.join(config['porechop_abi_path'], "{sample}_chop.fastq.gz")
    output:
        fasta=os.path.join(config['assemblies_work_path_shallow'], "{sample}","flye", "assembly.fasta")
    params:
        genome_size = lambda wildcards: genome_sizes[wildcards.sample],
        queue = "Priority",
        log = config['log'],
        outputdir = os.path.join(config['assemblies_work_path_shallow'], "{sample}","flye"),
        threads=12
    shell:
        """
        mkdir {params.outputdir} -p
        conda activate flye_2.9
        flye --nano-raw {input.fastq} --genome-size {params.genome_size} --out-dir {params.outputdir} --threads {params.threads} &&
        conda deactivate
        """
# Rule for assembling with Raven on shallow sequenced samples
rule raven_shallow:
    input:
        fastq=os.path.join(config['porechop_abi_path'], "{sample}_chop.fastq.gz")
    output:
        fasta=os.path.join(config['assemblies_work_path_shallow'], "{sample}","raven", "assembly.fasta")
    params:
        genome_size = lambda wildcards: genome_sizes[wildcards.sample],
        queue = "Priority",
        log = config['log'],
        outputdir = os.path.join(config['assemblies_work_path_shallow'], "{sample}","raven"),
        threads=12
    shell:
        """
        mkdir {params.outputdir} -p
        conda activate raven_env
        raven {input.fastq} --threads {params.threads} > {output.fasta}
        conda deactivate
        """

# Rule for assembling with Unicycler using subsamples
rule unicycler_sub:
    input:
        fastq=os.path.join(config['porechop_abi_path'], "{sample}_chop.fastq.gz")
    output:
        fasta=os.path.join(config['assemblies_work_path_shallow'], "{sample}","unicycler", "assembly.fasta")
    params:
        genome_size = lambda wildcards: genome_sizes[wildcards.sample],
        queue = "Priority",
        log = config['log'],
        outputdir = os.path.join(config['assemblies_work_path_shallow'], "{sample}","unicycler"),
        threads=12,
        keep = 0
    shell:
        """
        mkdir {params.outputdir} -p
        conda activate unicycler_0.5.0
        unicycler --long {input.fastq} --out {params.outputdir} --keep {params.keep} --threads {params.threads}
        conda deactivate
        """
