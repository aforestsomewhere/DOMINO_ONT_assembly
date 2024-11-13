# DOMINO WP3-4-Assemblies-ONT
# Author: Katie O'Mahony
# Adapted from: Julien Tap

shell.executable(config['shell_executable'])
shell.prefix(config['shell_prefix'])
genome_sizes = config['genome_sizes']

import os
import glob
import warnings
import gzip
#need color palette for nanocomp
import matplotlib
from matplotlib.colors import CSS4_COLORS

#Helper functions
# Function to generate a list of colors based on the number of barcodes
available_colours = list(CSS4_COLORS.keys())
def generate_colours(n):
    if n <= len(available_colours):
        return available_colours[:n]
    else:
        return (available_colours * ((n // len(available_colours)) + 1))[:n]

# Define the workflow
rule all:
    input:
        expand(os.path.join(config['porechop_abi_path'], "{sample}_chop.fastq.gz"), sample=config['ont_fastq_sample_names']),
        expand(os.path.join(config['filtlong_path'], "{sample}_filt.fastq.gz"), sample=config['ont_fastq_sample_names']),
        expand(os.path.join(config['nanoplot_path'], "{sample}", "NanoStats.txt"), sample=config['ont_fastq_sample_names']),
        os.path.join(config['nanocomp_path'], "NanoComp-report.html"),
        expand(os.path.join(config['flye_rough_path'], "{sample}", "flye01","assembly.fasta"),sample=config['ont_fastq_sample_names'])

# Rule for running Porechop_ABI on raw fastq files
rule porechop:
    input:
        os.path.join(config['ont_fastq_path'], "{sample}.fastq.gz")
    output:
        os.path.join(config['porechop_abi_path'], "{sample}_chop.fastq.gz")
    threads:
        config['porechop_threads']
    params:
        log=config['log'],
        queue="Priority"
    shell:
        "activate porechop_abi_env && porechop_abi --ab_initio -i {input} -o {output} --threads {threads} && conda deactivate"

# Rule for running Filtlong on fastq files post-porechop-abi
rule filtlong:
   input:
       os.path.join(config['porechop_abi_path'], "{sample}_chop.fastq.gz")
   output:
       os.path.join(config['filtlong_path'], "{sample}_filt.fastq.gz")
   threads:
       config['filtlong_threads']
   params:
       log = config['log'],
       queue = "Priority"
   shell:
       "activate filtlong_env && filtlong --min_length 1000 --keep_percent 95 {input} | gzip > {output} && conda deactivate"

# Rule for running NanoPlot on each fastq.gz file
# note - pandas warning: https://github.com/ranaroussi/yfinance/issues/1837
warnings.filterwarnings("ignore",
                        message="'T' is deprecated and will be removed in a future version, please use 'min' instead.",
                        category=FutureWarning, module="nanoplot")

rule nanoplot:
    input:
        os.path.join(config['filtlong_path'], "{sample}_filt.fastq.gz")
    output:
        os.path.join(config['nanoplot_path'], "{sample}", "NanoStats.txt"),
        os.path.join(config['nanoplot_path'], "{sample}", "NanoPlot-report.html")
    threads:
        config['nanoplot_threads']
    params:
        queue="Priority",
        log=config['log'],
        outdir=os.path.join(config['nanoplot_path'], "{sample}")
    shell:
        "activate nanocomp_env && NanoPlot --threads {threads} --fastq {input} --outdir {params.outdir} && conda deactivate"

# Rule for running NanoComp on the each fastq file
# note - pandas warning: https://github.com/ranaroussi/yfinance/issues/1837
warnings.filterwarnings("ignore",
                        message="'T' is deprecated and will be removed in a future version, please use 'min' instead.",
                        category=FutureWarning, module="nanocomp")

rule nanocomp:
    input:
	    expand(os.path.join(config['filtlong_path'],"{sample}_filt.fastq.gz"), sample=config['ont_fastq_sample_names'])
    output:
        os.path.join(config['nanocomp_path'], "NanoComp-report.html")
    params:
        queue="Priority",
        log=config['log'],
        outdir=config['nanocomp_path'],
        colours=generate_colours(len(config['ont_fastq_sample_names'])),
        names=[os.path.basename(f).replace('_filt.fastq.gz', '') for f in expand(os.path.join(config['filtlong_path'], "{sample}_filt.fastq.gz"), sample=config['ont_fastq_sample_names'])]
    shell:
        "activate nanocomp_env && NanoComp --threads {threads} --tsv_stats --outdir {params.outdir} --fastq {input} --names {params.names} --colors {params.colours} && conda deactivate"

rule flye_rough:
    input:
        fastq=os.path.join(config['filtlong_path'],"{sample}_filt.fastq.gz")
    output:
        fasta=os.path.join(config['flye_rough_path'], "{sample}","flye01", "assembly.fasta")
    params:
        genome_size = lambda wildcards: genome_sizes[wildcards.sample],
        queue = "Priority",
        log = config['log'],
        outputdir = os.path.join(config['flye_rough_path'], "{sample}","flye01"),
        threads=12
    shell:
        """
            mkdir {params.outputdir} -p && 
            conda activate flye_2.9 &&
            flye --nano-raw {input.fastq} --genome-size {params.genome_size} --out-dir {params.outputdir} --threads {params.threads} &&
            conda deactivate
        """
