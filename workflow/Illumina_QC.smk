# DOMINO WP3-4-Assemblies
# Author: Katie O'Mahony
# Adapted from: Julien Tap

shell.executable(config['shell_executable'])
shell.prefix(config['shell_prefix'])
genome_sizes = config['genome_sizes']

import os
import glob

# Define the workflow

#handle paired end reads
reads = ['1', '2']
rule all:
    input:
        config['multiqc_path'] + "/" + "multiqc_report.html",
        expand(os.path.join(config['fastqc_path'], "{sample}_R{read}_fastqc.html"), sample=config['illumina_fastq_sample_names'], read=reads),
        expand(os.path.join(config['mash_path'], "{sample}_R{read}.tab"), sample=config['illumina_fastq_sample_names'], read=reads),
        expand(os.path.join(config['fastp_path'], "{sample}_R{read}.fastq.gz"), sample=config['illumina_fastq_sample_names'], read=reads),
        expand(os.path.join(config['spades_path'], "{sample}_spades", "contigs.fasta"), sample=config['illumina_fastq_sample_names']),
        expand(os.path.join(config['spades_assemblies_path'], "{sample}_spades.fasta"), sample=config['illumina_fastq_sample_names'])

# Rule for running FastQC on raw fastq files
rule fastqc_raw:
    input:
        os.path.join(config['illumina_fastq_path'],"{sample}_R{read}.fastq.gz")
    output:
        os.path.join(config['fastqc_path'], "{sample}_R{read}_fastqc.html")
    params:
        queue = "Priority",
        log = config['log'],
        outputdirname = config['fastqc_path']
    shell:
        "activate fastqc_test && fastqc --threads 8 {input} -o {params.outputdirname} && conda deactivate"

# Rule for running MultiQC on FastQC reports
rule multiqc:
    input:
        expand(os.path.join(config['fastqc_path'], "{sample}_R{read}_fastqc.html"), sample=config['illumina_fastq_sample_names'], read=reads)
    output:
        report = config['multiqc_path'] + "/" + "multiqc_report.html"
    params:
        outputdir = config['multiqc_path'],
        inputdir = config['fastqc_path'],
        log = config['log'],
        queue = "Priority"
    shell:
        "activate multiqc_test && multiqc --outdir {params.outputdir} {params.inputdir} && conda deactivate"

# Rule for running Mash to screen readsets
rule mash:
    input:
        os.path.join(config['illumina_fastq_path'],"{sample}_R{read}.fastq.gz")
    output:
        os.path.join(config['mash_path'], "{sample}_R{read}.tab")
    params:
        outputdir = config['mash_path'],
        mash_db = config['mash_db'],
        log = config['log'],
        queue = "Priority"
    shell:
        "activate mash_test && mash screen -w -p 12 {params.mash_db} {input} > {output} && conda deactivate"

# Rule for running Fastp for quality filtering and trimming of screen readsets
rule fastp:
    input:
        r1 = os.path.join(config['illumina_fastq_path'],"{sample}_R1.fastq.gz"),
        r2 = os.path.join(config['illumina_fastq_path'],"{sample}_R2.fastq.gz")
    output:
        r1 = os.path.join(config['fastp_path'],"{sample}_R1.fastq.gz"),
        r2 = os.path.join(config['fastp_path'],"{sample}_R2.fastq.gz"),        
        json=os.path.join(config['fastp_path'],"{sample}_fastp.json"),
        html=os.path.join(config['fastp_path'],"{sample}_fastp.html")
    params:
        log = config['log'],
        queue = "Priority"
    shell:
        """
            conda activate fastp_env
            fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --json {output.json} --html {output.html}
            conda deactivate
        """

# Rule for assembling isolate reads using Spades
rule spades:
    input:
        r1 = os.path.join(config['fastp_path'],"{sample}_R1.fastq.gz"),
        r2 = os.path.join(config['fastp_path'],"{sample}_R2.fastq.gz")    
    output:
        assembly = os.path.join(config['spades_path'], "{sample}_spades", "contigs.fasta"),
        outdir = os.path.join(config['spades_path'], "{sample}_spades")
    params:
        queue = "Priority"
    shell:
        "activate SPAdes && spades.py -o {output.outdir} -1 {input.r1} -2 {input.r2} --only-assembler --isolate 2>&1 | tee -a {log} && conda deactivate"

#Rule to cleanup intermediate files
rule spades_cleanup:
    input:
        spadesfile = os.path.join(config['spades_path'], "{sample}_spades", "contigs.fasta")
    output:
        assembled = os.path.join(config['spades_assemblies_path'], "{sample}_spades.fasta")
    params:
        queue = "Priority"
    shell:
        """
        cp {input.spadesfile} {output.assembled}
        """
