# Main entrypoint of the workflow. 
# Main entrypoint of the workflow.
# Please follow the best practices:
# https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html,
# in particular regarding the standardized folder structure mentioned there.
shell.executable("/bin/bash")
shell.prefix("source /install/software/anaconda3.6.b/bin/activate;")

import os
import glob
import warnings

fastq_path = "/data/Food/analysis/R6564_NGS/Katie_F/DOMINO_WP4/test_ONT/fastq"

# Define the workflow
rule all:
    input:
        expand(os.path.join(config['nanoplot_path'],"{sample}", "NanoStats.txt"),sample=config['fastq_sample_names']),
        os.path.join(config['nanocomp_path'], "NanoStats.txt"),
        expand(os.path.join(config['porechop_path'], "{sample}_chop.fastq.gz"),sample=config['fastq_sample_names'])

# Rule for running NanoPlot on each fastq.gz file
#note - pandas warning: https://github.com/ranaroussi/yfinance/issues/1837
warnings.filterwarnings("ignore", message="'T' is deprecated and will be removed in a future version, please use 'min' instead.", category=FutureWarning, module="nanoplot")

rule nanoplot:
    input:
        fastq_path + "/" + "{sample}.fastq.gz"
    output:
        os.path.join(config['nanoplot_path'],"{sample}","NanoStats.txt"),
        os.path.join(config['nanoplot_path'],"{sample}","NanoPlot-report.html")
    threads:
        config['nanoplot_threads']
    params:
        queue = "Priority",
        log = config['log'],
        outdir = os.path.join(config['nanoplot_path'],"{sample}")
    shell:
        "conda activate nanocomp_env && NanoPlot --threads {threads} --fastq {input} --outdir {params.outdir} && conda deactivate"

# Rule for running NanoComp on overall sequencing summary file
#note - pandas warning: https://github.com/ranaroussi/yfinance/issues/1837
warnings.filterwarnings("ignore", message="'T' is deprecated and will be removed in a future version, please use 'min' instead.", category=FutureWarning, module="nanocomp")

rule nanocomp:
    input:
        config['sequencing_summary']
    output:
        os.path.join(config['nanocomp_path'],"NanoStats.txt"),
        os.path.join(config['nanocomp_path'],"NanoComp-report.html")
    threads:
        config['nanocomp_threads']
    params:
        queue = "Priority",
        log = config['log'],
        outdir = config['nanocomp_path']
    shell:
        "conda activate nanocomp_env && NanoComp --threads {threads} --tsv_stats --outdir {params.outdir} --barcoded --summary {input} && conda deactivate"

# Rule for running Porechop_ABI on raw fastq files
rule porechop:
   input:
       fastq_path + "/" + "{sample}.fastq.gz"
   output:
       os.path.join(config['porechop_path'], "{sample}_chop.fastq.gz")
   threads:
       config['porechop_threads']
   params:
       log = config['log'],
       queue = "Priority"
   shell:
       "conda activate porechop_abi_env && porechop_abi --ab_initio -i {input} -o {output} --threads {threads} && conda deactivate"
