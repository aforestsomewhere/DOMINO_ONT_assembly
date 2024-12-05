# DOMINO WP3-4-Assemblies-ONT
# Assembly for adequate-depth isolates, with sub-sampling
# Author: Katie O'Mahony
# Adapted from: Julien Tap

shell.executable(config['shell_executable'])
shell.prefix(config['shell_prefix'])

import os
import glob
import warnings
import gzip

genome_sizes = config['genome_sizes']
###HELPER FUNCTIONS #####

# Define the workflow
rule all:
    input:
        expand(os.path.join(config['subsamples_path'], "{sample}", "sample_0{num}.fastq"),sample=config['deep_sample_names'],num=range(1, 4)),
         #subsampled assemblies
        expand(os.path.join(config['assemblies_work_path'], "{sample}", "flye0{num}","assembly.fasta"),sample=config['deep_sample_names'],num=range(1, 4)),
        expand(os.path.join(config['assemblies_work_path'], "{sample}", "raven0{num}","assembly.fasta"),sample=config['deep_sample_names'],num=range(1, 4)),
        expand(os.path.join(config['assemblies_work_path'], "{sample}", "unicycler0{num}","assembly.fasta"),sample=config['deep_sample_names'],num=range(1, 4))

# Rule for subsampling reads using Trycycler subsample
rule subsampling:
    input:
        os.path.join(config['filtlong_path'],"{sample}_filt.fastq.gz")
    output:
        os.path.join(config['subsamples_path'], "{sample}", "sample_0{num}.fastq")
    params:
        genome_size = lambda wildcards: genome_sizes[wildcards.sample],
        log = config['log'],
        outputdir = directory(os.path.join(config['subsamples_path'], "{sample}"))
    threads: 8
    shell:
        """
        conda activate trycycler_env
        trycycler subsample --reads {input} --out_dir {params.outputdir} --count 3 --genome_size {params.genome_size} --min_read_depth 50 --threads {threads}
        conda deactivate
        """

# Rule for assembling with flye using subsamples
rule flye_sub:
    input:
        fastq=os.path.join(config['subsamples_path'], "{sample}", "sample_0{num}.fastq"),
        shallow=os.path.join(config['filtlong_path'],"{sample}_filt.fastq.gz")
    output:
        fasta=os.path.join(config['assemblies_work_path'], "{sample}","flye0{num}", "assembly.fasta")
    params:
        genome_size = lambda wildcards: genome_sizes[wildcards.sample],
        queue = "Priority",
        log = config['log'],
        outputdir = os.path.join(config['assemblies_work_path'], "{sample}","flye0{num}"),
        outputdirshallow = os.path.join(config['assemblies_work_path'], "{sample}","flye01"),
        threads=12
    shell:
        """
        if [[ -f {input.fastq} ]]; then
            mkdir {params.outputdir} -p && 
            echo {input.fastq} >> subsample_ok.txt &&
            conda activate flye_2.9 &&
            flye --nano-raw {input.fastq} --genome-size {params.genome_size} --out-dir {params.outputdir} --plasmids --threads {params.threads} &&
            conda deactivate
        else
            mkdir {params.outputdir} -p && 
            conda activate flye_2.9 &&
            flye --nano-raw {input.shallow} --genome-size {params.genome_size} --out-dir {params.outputdirshallow} --plasmids --threads {params.threads} &&
            conda deactivate
        fi
        """

# Rule for assembling with Raven using subsamples
rule raven_sub:
    input:
        fastq=os.path.join(config['subsamples_path'], "{sample}", "sample_0{num}.fastq")
    output:
        fasta=os.path.join(config['assemblies_work_path'], "{sample}","raven0{num}", "assembly.fasta")
    params:
        genome_size = lambda wildcards: genome_sizes[wildcards.sample],
        queue = "Priority",
        log = config['log'],
        outputdir = os.path.join(config['assemblies_work_path'], "{sample}","raven0{num}"),
        threads=12
    shell:
        """
        if [[ -f {input.fastq} ]]; then
            mkdir {params.outputdir} -p
            conda activate raven_env && \
            raven {input.fastq} --threads {params.threads} > {output.fasta} && \
            conda deactivate
        else
            touch {output}
        fi
        """

# Rule for assembling with Unicycler using subsamples
rule unicycler_sub:
    input:
        fastq=os.path.join(config['subsamples_path'], "{sample}", "sample_0{num}.fastq")
    output:
        fasta=os.path.join(config['assemblies_work_path'], "{sample}","unicycler0{num}", "assembly.fasta")
    params:
        genome_size = lambda wildcards: genome_sizes[wildcards.sample],
        queue = "Priority",
        log = config['log'],
        outputdir = os.path.join(config['assemblies_work_path'], "{sample}","unicycler0{num}"),
        threads=12,
        keep = 0
    shell:
        """
        if [[ -f {input.fastq} ]]; then
            mkdir {params.outputdir} -p
            conda activate unicycler_0.5.0
            unicycler --long {input.fastq} --out {params.outputdir} --keep {params.keep} --threads {params.threads} && conda deactivate
        else
            touch {output}
        fi
        """















