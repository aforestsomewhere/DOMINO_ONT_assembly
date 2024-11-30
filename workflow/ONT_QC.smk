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
        expand(os.path.join(config['mash_path'], "{sample}.tab"), sample=config['ont_fastq_sample_names']),
        expand(os.path.join(config['mash_path'], "{sample}.png"), sample=config['ont_fastq_sample_names']),
        expand(os.path.join(config['porechop_abi_path'], "{sample}_chop.fastq.gz"), sample=config['ont_fastq_sample_names']),
        expand(os.path.join(config['filtlong_path'], "{sample}_filt.fastq.gz"), sample=config['ont_fastq_sample_names']),
        expand(os.path.join(config['nanoplot_path'], "{sample}", "NanoStats.txt"), sample=config['ont_fastq_sample_names']),
        os.path.join(config['nanocomp_path'], "NanoComp-report.html"),
        expand(os.path.join(config['flye_rough_path'], "{sample}", "flye01","assembly.fasta"),sample=config['ont_fastq_sample_names'])

# Rule for profiling readsets against Mash DB to detect potential contamination/leakage
rule mash:
    input:
        os.path.join(config['ont_fastq_path'], "{sample}.fastq.gz")
    output:
        os.path.join(config['mash_path'], "{sample}.tab")
    log:
        out =os.path.join(config['log'],"{sample}_mash_bac.out") ,
        err = os.path.join(config['log'],"{sample}_mash_bac.err")
    params:
        mash_db = config['mash_db'],
        log=config['log'],
        queue = "Priority"
    shell:
        "activate mash_test && mash screen -w -p 12 {params.mash_db} {input} 1> {log.out} 2> {log.err} > {output} && conda deactivate"

rule mash_plots:
    input:
        os.path.join(config['mash_path'], "{sample}.tab")
    output:
        os.path.join(config['mash_path'], "{sample}.png")
    params:
        width=80,  # Specify wrap width as a parameter for flexibility
        colors=[
            "#432a74", "#8c61fc", "#00babe", "#ff8c00", "#c7c4fa", "#bcebeb", "#ffe6cc", "#ff7f00", "#984ea3", "#ffff33"
        ]  # Custom colors
    run:
        import pandas as pd
        import matplotlib.pyplot as plt

        # Define a function to wrap text by inserting newline characters
        def wrap_text(text, width):
            return '\n'.join([text[i:i+width] for i in range(0, len(text), width)])

        # Load data
        df = pd.read_csv(input[0], sep='\t', header=None)

        # Rename columns
        df.columns = ["identity", "shared-hashes", "median-multiplicity", "p-value", "query-ID", "query-description"]

        # Extract match counts
        df['match_count'] = df['shared-hashes'].str.split('/').str[0].astype(int)

        # Get top entries based on match count
        top_df = df.sort_values(by="match_count", ascending=False).head(10)

        # Apply wrapping to the 'query-description' column for y-axis labels
        top_df['wrapped_query-description'] = top_df['query-description'].apply(lambda x: wrap_text(x, width=params.width))

        # Plot the horizontal bar chart
        plt.figure(figsize=(14, 8))
        plt.barh(top_df['wrapped_query-description'], top_df['match_count'], color=params.colors)
        plt.xlabel('Shared K-mer Matches')
        plt.ylabel('Genomic Reference')
        plt.title('Top MASH Screen Shared K-mer Matches')

        # Adjust layout to make room for y-axis labels
        plt.subplots_adjust(left=0.6)
        plt.gca().invert_yaxis()

        # Save the plot to the output path
        plt.savefig(output[0], format='png', dpi=300)
        plt.close()

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
        fasta=os.path.join(config['flye_rough_path'], "{sample}","flye01", "assembly.fasta"),
        logfile=os.path.join(config['flye_rough_path'], "{sample}", "flye01", "flye.log")
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
