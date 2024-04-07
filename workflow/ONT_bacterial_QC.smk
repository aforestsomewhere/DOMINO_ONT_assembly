# DOMINO WP3-4-Assemblies-ONT
# Author: Katie O'Mahony
# Adapted from: Julien Tap

shell.executable(config['shell_executable'])
shell.prefix(config['shell_prefix'])

import os
import glob
import warnings
import gzip

# fastq_path = config['raw_fastq_path']

###HELPER FUNCTIONS #####

# Define the workflow
rule all:
    input:
        expand(os.path.join(config['nanoplot_path'], "{sample}", "NanoStats.txt"), sample=config['fastq_sample_names']),
        os.path.join(config['nanocomp_path'], "NanoStats.txt"),
        expand(os.path.join(config['porechop_path'], "{sample}_chop.fastq.gz"), sample=config['fastq_sample_names']),
        expand(os.path.join(config['filtlong_path'], "{sample}_filt.fastq.gz"), sample=config['fastq_sample_names'])

# Rule for running NanoPlot on each fastq.gz file
# note - pandas warning: https://github.com/ranaroussi/yfinance/issues/1837
warnings.filterwarnings("ignore",
                        message="'T' is deprecated and will be removed in a future version, please use 'min' instead.",
                        category=FutureWarning, module="nanoplot")

rule nanoplot:
    input:
        os.path.join(config['raw_fastq_path'], "{sample}.fastq.gz")
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
        "conda activate nanocomp_env && NanoPlot --threads {threads} --fastq {input} --outdir {params.outdir} && conda deactivate"

# Rule for running NanoComp on the overall sequencing summary file
# note - pandas warning: https://github.com/ranaroussi/yfinance/issues/1837
warnings.filterwarnings("ignore",
                        message="'T' is deprecated and will be removed in a future version, please use 'min' instead.",
                        category=FutureWarning, module="nanocomp")

rule nanocomp:
    input:
        config['sequencing_summary']
    output:
        os.path.join(config['nanocomp_path'], "NanoStats.txt"),
        os.path.join(config['nanocomp_path'], "NanoComp-report.html")
    threads:
        config['nanocomp_threads']
    params:
        queue="Priority",
        log=config['log'],
        outdir=config['nanocomp_path']
    shell:
        "conda activate nanocomp_env && NanoComp --threads {threads} --tsv_stats --outdir {params.outdir} --barcoded --summary {input} && conda deactivate"

# Rule for running Porechop_ABI on raw fastq files
rule porechop:
    input:
        os.path.join(config['raw_fastq_path'], "{sample}.fastq.gz")
    output:
        os.path.join(config['porechop_path'], "{sample}_chop.fastq.gz")
    threads:
        config['porechop_threads']
    params:
        log=config['log'],
        queue="Priority"
    shell:
        "conda activate porechop_abi_env && porechop_abi --ab_initio -i {input} -o {output} --threads {threads} && conda deactivate"

# Rule for running Filtlong on fastq files post-porechop-abi
rule filtlong:
   input:
       os.path.join(config['porechop_path'], "{sample}_chop.fastq.gz")
   output:
       os.path.join(config['filtlong_path'], "{sample}_filt.fastq.gz")
   threads:
       config['filtlong_threads']
   params:
       log = config['log'],
       queue = "Priority"
   shell:
       "conda activate filtlong_env && filtlong --min_length 1000 --keep_percent 95 {input} | gzip > {output} && conda deactivate"
