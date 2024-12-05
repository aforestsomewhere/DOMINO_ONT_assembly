# DOMINO WP3-4-Assemblies-ONT
# Author: Katie O'Mahony
# Adapted from: Julien Tap

shell.executable(config['shell_executable'])
shell.prefix(config['shell_prefix'])

import os
import glob
import warnings
import gzip

rule all:
    input:
        expand(os.path.join(config['trycycler_path'],"{sample}","assemblies","{sample}_{assembler}.fasta"), sample=config['deep_sample_names'],assembler=config['assemblers']),
        expand(os.path.join(config['trycycler_path'],"{sample}","{sample}.fastq.gz"), sample=config['deep_sample_names']),
        expand(os.path.join(config['trycycler_path'], "{sample}", "trycycler", "contigs.phylip"),sample=config['deep_sample_names'])
       #  expand(os.path.join(config['trycycler_path'], "{sample}","clustering_ok.txt"), sample=config['subsampled_sample_names'])

rule tidy_assemblies:
    input:
        assembly = os.path.join(config['assemblies_work_path'], "{sample}", "{assembler}", "assembly.fasta")
    output:
        assembly = os.path.join(config['trycycler_path'],"{sample}","assemblies","{sample}_{assembler}.fasta")
    params:
        queue = "Priority",
        outdir = os.path.join(config['trycycler_path'], "{sample}","assemblies"),
        tryc = os.path.join(config['trycycler_path'])
    shell:
        """
            mkdir -p {params.tryc} && mkdir -p {params.outdir} && cp {input.assembly} {output.assembly}
            exit 0
        """
rule tidy_fastq:
    input:
        fastq = os.path.join(config['ont_fastq_path'], "{sample}.fastq.gz")
    output:
        fastq = os.path.join(config['trycycler_path'],"{sample}","{sample}.fastq.gz")
    params:
        queue = "Priority"
    shell:
        "cp {input.fastq} {output.fastq}"

rule cluster_subsamples:
    input:
        #assemblies = os.path.join(config['trycycler_path'], "{sample}", "assemblies"),
        assemblies = expand(os.path.join(config['trycycler_path'],"{sample}","assemblies","{sample}_{assembler}.fasta"),sample=config['deep_sample_names'],assembler=config['assemblers']),
        fastq = os.path.join(config['trycycler_path'],"{sample}", "{sample}.fastq.gz")
    output:
        os.path.join(config['trycycler_path'], "{sample}", "trycycler", "contigs.phylip")
        #os.path.join(config['trycycler_path'], "{sample}", "clustering_ok.txt")
    params:
        cluster_dir=os.path.join(config['trycycler_path'],"{sample}","trycycler"),
        threads=16,
        log = config['log'],
        queue = "Priority"
    shell:
      """
      rm -rf {params.cluster_dir}
      conda activate trycycler_env
      trycycler cluster --assemblies {input.assemblies} --threads {params.threads} --reads {input.fastq} --out_dir {params.cluster_dir}
      conda deactivate
      """
