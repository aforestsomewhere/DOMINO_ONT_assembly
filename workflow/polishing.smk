# DOMINO WP3-4-Assemblies-ONT
# Author: Katie O'Mahony
# Adapted from: Julien Tap
# Workflow - short read polishing
shell.executable(config['shell_executable'])
shell.prefix(config['shell_prefix'])

import os
import glob
import warnings
import gzip

# Define the workflow
rule all:
    input:
        fasta=expand(os.path.join(config['polished_assemblies'], "{sample}_polished.fasta"), sample=config['polishing_sample_names']),
        index = expand(os.path.join(config['assemblies_work_path'], "{sample}.fasta.ann"), sample=config['polishing_sample_names']),
        align_r1 = expand(os.path.join(config['assemblies_work_path'], "{sample}_alignments_1.sam"), sample=config['polishing_sample_names']),
        align_r2 = expand(os.path.join(config['assemblies_work_path'], "{sample}_alignments_2.sam"), sample=config['polishing_sample_names']),
        filt_r1 = expand(os.path.join(config['assemblies_work_path'], "{sample}_filtered_1.sam"), sample=config['polishing_sample_names']),
        filt_r2 = expand(os.path.join(config['assemblies_work_path'], "{sample}_filtered_2.sam"), sample=config['polishing_sample_names']),
        flag = expand(os.path.join(config['polished_assemblies'], "{sample}.flag"), sample=config['polishing_sample_names'])

coverage_values = config['genome_coverage']
#Rules for any samples which had insufficent depth for subsampling
# Rule for assembling with flye without subsamples
#use same environment for pypolca and polypolish (which needs bwa, pypolca_env has version 0.7.18)
#working on command line - running out of mem when run together? 
rule bwa_index:
    input:
        fasta = os.path.join(config['assemblies_work_path'], "{sample}.fasta")
    output:
        index = os.path.join(config['assemblies_work_path'], "{sample}.fasta.ann")
    params:
        coverage = lambda wildcards: coverage_values[wildcards.sample],
        log = config['log'],
        threads=config['pypolca_threads']
    shell:
        "activate pypolca_env && bwa index {input.fasta} && conda deactivate"

rule bwa_mem_1:
    input:
        fasta = os.path.join(config['assemblies_work_path'], "{sample}.fasta"),
        r1 = os.path.join(config['illumina_fastq_path'], "{sample}_R1.fastq.gz")
    output:
        align_r1 = os.path.join(config['assemblies_work_path'], "{sample}_alignments_1.sam")
    params:
        log = config['log'],
        outputdir = os.path.join(config['assemblies_work_path'], "{sample}","trycycler_pypolca"),
        threads=config['pypolca_threads']
    shell:
        "activate pypolca_env && bwa mem -t {params.threads} -a {input.fasta} {input.r1} > {output.align_r1} && conda deactivate"
rule bwa_mem_2:
    input:
        fasta = os.path.join(config['assemblies_work_path'], "{sample}.fasta"),
        r2 = os.path.join(config['illumina_fastq_path'], "{sample}_R2.fastq.gz")
    output:
        align_r2 = os.path.join(config['assemblies_work_path'], "{sample}_alignments_2.sam")
    params:
        coverage = lambda wildcards: coverage_values[wildcards.sample],
        log = config['log'],
        outputdir = os.path.join(config['assemblies_work_path'], "{sample}","trycycler_pypolca"),
        threads=config['pypolca_threads']
    shell:
        "activate pypolca_env && bwa mem -t {params.threads} -a {input.fasta} {input.r2} > {output.align_r2} && conda deactivate"
rule polypolish_filter:
    input:
        align_r1 = os.path.join(config['assemblies_work_path'], "{sample}_alignments_1.sam"),
        align_r2 = os.path.join(config['assemblies_work_path'], "{sample}_alignments_2.sam")
    output:
        filt_r1 = os.path.join(config['assemblies_work_path'], "{sample}_filtered_1.sam"),
        filt_r2 = os.path.join(config['assemblies_work_path'], "{sample}_filtered_2.sam")
    params:
        log = config['log'],
        threads=config['pypolca_threads'],
        polypolish = config['polypolish_binary']
    shell:
        "activate && {params.polypolish} filter --in1 {input.align_r1} --in2 {input.align_r2} --out1 {output.filt_r1} --out2 {output.filt_r2}"

rule polypolish_polish:
    input:
        fasta = os.path.join(config['assemblies_work_path'], "{sample}.fasta"),
        r1 = os.path.join(config['illumina_fastq_path'], "{sample}_R1.fastq.gz"),       
        r2 = os.path.join(config['illumina_fastq_path'], "{sample}_R2.fastq.gz"),
        filt_r1 = os.path.join(config['assemblies_work_path'], "{sample}_filtered_1.sam"),
        filt_r2 = os.path.join(config['assemblies_work_path'], "{sample}_filtered_2.sam")
    output:
        fasta=os.path.join(config['polished_assemblies'], "{sample}_polished.fasta"),
        pypolcareport = os.path.join(config['polished_assemblies'], "{sample}.report"),
        flag = os.path.join(config['polished_assemblies'], "{sample}.flag")
    params:
        coverage = lambda wildcards: coverage_values[wildcards.sample],
        name = lambda wildcards: wildcards.sample,
        tempfasta = os.path.join(config['assemblies_work_path'],"{sample}_poly.fasta"),
        polypolish = config['polypolish_binary'],
        log = config['log'],
        outputdir = os.path.join(config['assemblies_work_path'], "{sample}","trycycler_pypolca"),
        threads=config['pypolca_threads']
    shell:
        """
        #convert to integer - bash doesn't like floating decimals here
        set -e
        set -x
        coverage_value=$(echo "{params.coverage}" | bc -l)
        coverage_value_rounded=$(echo "($coverage_value + 0.5) / 1" | bc)
        echo $coverage_value_rounded
        set +e
        set +x
        #Coverage <5X - only Polypolish careful
        if [[ $coverage_value_rounded -lt 5 ]]; then
            {params.polypolish} polish --careful {input.fasta} {input.filt_r1} {input.filt_r2} > {params.tempfasta}
            conda deactivate
            mv {params.tempfasta} {output.fasta}
            touch {output.flag}
        #Coverage <25X and >5X -  Polypolish careful + Pypolca careful
        elif [[ $coverage_value_rounded -ge 5 && $coverage_value_rounded -lt 25 ]]; then
            {params.polypolish} polish --careful {input.fasta} {input.filt_r1} {input.filt_r2} > {params.tempfasta}
            #also pypolca careful
            conda activate pypolca_env
            pypolca run -a {params.tempfasta} -1 {input.r1} -2 {input.r2} -p {params.name} --out-dir {params.outputdir} --threads {params.threads} --careful
            conda deactivate
            touch {output.flag}
        #Coverage >25X Polypolish (regular) + Pypolca careful
        elif [[ $coverage_value_rounded -gt 25 ]]; then
            {params.polypolish} polish {input.fasta} {input.filt_r1} {input.filt_r2} > {params.tempfasta}
            #also pypolca careful
            conda activate pypolca_env
            pypolca run -a {params.tempfasta} -1 {input.r1} -2 {input.r2} -p {params.name} --threads {params.threads} --careful -o {params.outputdir}
            cp {params.outputdir}/{params.name}_corrected.fasta {output.fasta}
            cp {params.outputdir}/{params.name}.report {output.pypolcareport}
            conda deactivate
            touch {output.flag}
        fi
        rm -rf {params.outputdir}
        """
