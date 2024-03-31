# DOMINO WP3-4-Assemblies-ONT
# Author: Katie O'Mahony
# Adapted from: Julien Tap

shell.executable("/bin/bash")
shell.prefix("source /install/software/anaconda3.6.b/bin/activate;")

import os
import glob
import warnings
import gzip

fastq_path = "/data/Food/analysis/R6564_NGS/Katie_F/DOMINO_WP4/test_ONT/fastq"

###HELPER FUNCTIONS #####

# Define a function to retrieve genome_size based on the sample name
def get_genome_size(wildcards):
    return config['genome_sizes'].get(wildcards.sample, None)

# Define a function to count the number of reads in a fastq.gz
def count_fastq_reads(file_path):
    with gzip.open(file_path, 'rt') as f:
        num_reads = sum(1 for _ in f) // 4
    return num_reads

#Define a function to count the number of bases in a fastq.gz
def count_fastq_bp(file_path):
    total_bp = 0
    with gzip.open(file_path, 'rt') as f:
        for line in f:
            # Exclude header lines starting with ">"
            if not line.startswith('>'):
                total_bp += len(line.strip())
    return total_bp

#Defines a function to generate updated sample list for thos with adequate depth for subsampling
def generate_input_list():
    input_list = []
    for sample in config['fastq_sample_names']:
        depth_path = os.path.join(config['depth_path'], f"{sample}.txt")
        with open(depth_path) as depth_file:
            depth_value = float(depth_file.read().strip())
        if depth_value > 25:
            input_list.append(sample)
    return input_list

# Define the workflow
rule all:
    input:
        expand(os.path.join(config['nanoplot_path'], "{sample}", "NanoStats.txt"), sample=config['fastq_sample_names']),
        os.path.join(config['nanocomp_path'], "NanoStats.txt"),
        expand(os.path.join(config['porechop_path'], "{sample}_chop.fastq.gz"), sample=config['fastq_sample_names']),
        expand(os.path.join(config['filtlong_path'], "{sample}_filt.fastq.gz"), sample=config['fastq_sample_names']),
        expand(os.path.join(config['depth_path'], "{sample}.txt"), sample=config['fastq_sample_names']),
        os.path.join(config['depth_path'], "inputlist.txt"),
        #subsamples
        expand(os.path.join(config['trycycler_subsample_path'], "{sample}", "sample_0{num}.fastq"),sample=config['fastq_sample_names'],num=range(1, 4)),
        #subsampled assemblies
        expand(os.path.join(config['assemblies_work_path'], "{sample}", "flye0{num}","assembly.fasta"),sample=config['fastq_sample_names'],num=range(1, 4)),
        expand(os.path.join(config['assemblies_work_path'], "{sample}", "raven0{num}","assembly.fasta"),sample=config['fastq_sample_names'],num=range(1, 4)),
        expand(os.path.join(config['assemblies_work_path'], "{sample}", "unicycler0{num}","assembly.fasta"),sample=config['fastq_sample_names'],num=range(1, 4)),
        #full assemblies (where no subsampling possible)
        expand(os.path.join(config['assemblies_work_path'], "{sample}", "flyefull","assembly.fasta"),sample=config['fastq_sample_names']),
        expand(os.path.join(config['assemblies_work_path'], "{sample}", "ravenfull","assembly.fasta"),sample=config['fastq_sample_names']),
        expand(os.path.join(config['assemblies_work_path'], "{sample}", "unicyclerfull","assembly.fasta"),sample=config['fastq_sample_names']),
        expand(os.path.join(config['assemblies_path'],"{sample}"),sample=config['fastq_sample_names']),
        expand(os.path.join(config['assemblies_path'],"{sample}","{sample}_{assembler}.fasta"), sample=config['fastq_sample_names'],assembler=config['assemblers']),
        expand(os.path.join(config['trycycler_cluster_path'], "{sample}","clustering_full_ok.txt"), sample=config['fastq_sample_names'])

#Rule for running NanoPlot on each fastq.gz file
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

# Rule for subsampling reads with Trycycler (post-porechop-abi, post-filtlong)
#only subsample at adequate depth (>25x), otherwise assemble directly from all reads
rule estimate_depth:
   input:
       os.path.join(config['filtlong_path'], "{sample}_filt.fastq.gz")
   output:
       os.path.join(config['depth_path'], "{sample}.txt")
   threads: 1
   params:
       log = config['log'],
       queue = "Priority",
       genome_size = get_genome_size 
#       depth = lambda wildcards: count_fastq_reads(input[0])
   run:
        total_bp = count_fastq_bp(input[0])
        genome_size_numeric = float(params.genome_size[:-1]) * (10 ** 6)  # Convert genome size from "X.Xm" to X.X
        depth = total_bp/genome_size_numeric
        with open(output[0], 'w') as out_file:
            out_file.write(str(depth))

rule generate_input_list:
    # Rule to generate input list
    input:
        # Depth files as input
        depth_files = expand(os.path.join(config['depth_path'], "{sample}.txt"), sample=config['fastq_sample_names'])
    output:
        os.path.join(config['depth_path'],"inputlist.txt")
    params:
        log = config['log'],
        queue = "Priority"
    run:
        input_list = generate_input_list()
        shell("""
                echo {input_list} > {output}
                """)
        # Store the generated input list to a file or use it directly in subsequent rules conda

rule subsample:
    input:
        fastq = os.path.join(config['filtlong_path'], "{sample}_filt.fastq.gz"),
        depth = os.path.join(config['depth_path'], "{sample}.txt")
    output:
        os.path.join(config['trycycler_subsample_path'], "{sample}", "sample_0{num}.fastq")
    threads:
        config['trycycler_subsample_threads']
    params:
        log = config['log'],
        genome_size = get_genome_size,
        queue = "Priority",
        outputdir = os.path.join(config['trycycler_subsample_path'], "{sample}")
    shell:
        """
        ARG=$(cat {input.depth})
        if (( $(echo "$ARG > 20" | bc -l) )); then
            conda activate trycycler_env &&
            mkdir -p {params.outputdir} &&
            trycycler subsample --reads {input.fastq} --out_dir {params.outputdir} --count 3 --genome_size {params.genome_size} --min_read_depth 20 &&
            conda deactivate
        else
            touch {output}
        fi
        """

# Rule for assembling with flye using subsamples
rule flye_sub:
    input:
        fastq=os.path.join(config['trycycler_subsample_path'], "{sample}", "sample_0{num}.fastq"),
        depth = os.path.join(config['depth_path'], "{sample}.txt")
    output:
        fasta=os.path.join(config['assemblies_work_path'], "{sample}","flye0{num}", "assembly.fasta")
    params:
        genome_size = get_genome_size,
        queue = "Priority",
        log = config['log'],
        outputdir = os.path.join(config['assemblies_work_path'], "{sample}","flye0{num}"),
        threads=config['flye_threads']
    shell:
        """
        ARG=$(cat {input.depth})
        if (( $(echo "$ARG > 20" | bc -l) )); then
            mkdir {params.outputdir} -p && 
            conda activate flye_2.9 &&
            flye --nano-raw {input.fastq} --genome-size {params.genome_size} --out-dir {params.outputdir} --plasmids --threads {params.threads} &&
            conda deactivate
        else
            touch {output}
        fi
        """

# Rule for assembling with Raven using subsamples
rule raven_sub:
    input:
        fastq=os.path.join(config['trycycler_subsample_path'], "{sample}", "sample_0{num}.fastq"),
        depth = os.path.join(config['depth_path'], "{sample}.txt")
    output:
        fasta=os.path.join(config['assemblies_work_path'], "{sample}","raven0{num}", "assembly.fasta")
    params:
        genome_size = get_genome_size,
        queue = "Priority",
        log = config['log'],
        outputdir = os.path.join(config['assemblies_work_path'], "{sample}","raven0{num}"),
        threads=config['raven_threads']
    shell:
        """
        ARG=$(cat {input.depth})
        if (( $(echo "$ARG > 20" | bc -l) )); then
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
        fastq=os.path.join(config['trycycler_subsample_path'], "{sample}", "sample_0{num}.fastq"),
        depth = os.path.join(config['depth_path'], "{sample}.txt")
    output:
        fasta=os.path.join(config['assemblies_work_path'], "{sample}","unicycler0{num}", "assembly.fasta")
    params:
        genome_size = get_genome_size,
        queue = "Priority",
        log = config['log'],
        outputdir = os.path.join(config['assemblies_work_path'], "{sample}","unicycler0{num}"),
        threads=config['unicycler_threads'],
        keep = 0
    shell:
        """
        ARG=$(cat {input.depth})
        if (( $(echo "$ARG > 20" | bc -l) )); then
            mkdir {params.outputdir} -p
            conda activate unicycler_0.5.0 && \
            unicycler --long {input.fastq} --out {params.outputdir} --keep {params.keep} --threads {params.threads} && conda deactivate
        else
            touch {output}
        fi
        """


#Rules for any samples which had insufficent depth for subsampling
# Rule for assembling with flye without subsamples
rule flye_full:
    input:
        fastq = os.path.join(config['filtlong_path'], "{sample}_filt.fastq.gz"),
        depth = os.path.join(config['depth_path'], "{sample}.txt")
    output:
        fasta=os.path.join(config['assemblies_work_path'], "{sample}","flyefull", "assembly.fasta")
    params:
        genome_size = get_genome_size,
        queue = "Priority",
        log = config['log'],
        outputdir = os.path.join(config['assemblies_work_path'], "{sample}","flyefull"),
        threads=config['flye_threads']
    shell:
        """
        mkdir {params.outputdir} -p && 
        conda activate flye_2.9 &&
        flye --nano-raw {input.fastq} --genome-size {params.genome_size} --out-dir {params.outputdir} --plasmids --threads {params.threads} &&
        conda deactivate
        touch {output}
        """

# Rule for assembling with raven without subsamples
rule raven_full:
    input:
        fastq = os.path.join(config['filtlong_path'], "{sample}_filt.fastq.gz"),
        depth = os.path.join(config['depth_path'], "{sample}.txt")
    output:
        fasta=os.path.join(config['assemblies_work_path'], "{sample}","ravenfull", "assembly.fasta")
    params:
        genome_size = get_genome_size,
        queue = "Priority",
        log = config['log'],
        outputdir = os.path.join(config['assemblies_work_path'], "{sample}","ravenfull"),
        threads=config['raven_threads']
    shell:
        """
        mkdir {params.outputdir} -p && \
        conda activate raven_env && \
        raven {input.fastq} --threads {params.threads} > {output.fasta} && \
        conda deactivate
        touch {output}
        """

# Rule for assembling with Unicycler without subsamples
rule unicycler_full:
    input:
        fastq = os.path.join(config['filtlong_path'], "{sample}_filt.fastq.gz"),
        depth = os.path.join(config['depth_path'], "{sample}.txt")
    output:
        fasta=os.path.join(config['assemblies_work_path'], "{sample}","unicyclerfull", "assembly.fasta")
    params:
        genome_size = get_genome_size,
        queue = "Priority",
        log = config['log'],
        outputdir = os.path.join(config['assemblies_work_path'], "{sample}","unicyclerfull"),
        threads=config['unicycler_threads'],
        keep = 0
    shell:
        """
        mkdir {params.outputdir} -p
        conda activate unicycler_0.5.0 &&
        unicycler --long {input.fastq} --out {params.outputdir} --keep {params.keep} --threads {params.threads} || true && conda deactivate
        if [[ ! -f {output} ]]; then
            touch {output}
        fi
        """
rule tidy_assemblies:
    input:
        assembly = os.path.join(config['assemblies_work_path'], "{sample}", "{assembler}", "assembly.fasta")
    output:
        assembly = os.path.join(config['assemblies_path'], "{sample}", "{sample}_{assembler}.fasta")
    params:
        queue = "Priority",
        outdir = os.path.join(config['assemblies_path'], "{sample}")
    shell:
        """
        mkdir -p {params.outdir} && cp {input.assembly} {output.assembly}
        exit 0
        """
# Rule for Trycycler clustering
rule clustering_full:
    input:
        assembly=os.path.join(config['assemblies_path'], "{sample}"),
        fastq=os.path.join(config['filtlong_path'], "{sample}_filt.fastq.gz"),
        depth=os.path.join(config['depth_path'], "{sample}.txt")
    output:
        os.path.join(config['trycycler_cluster_path'], "{sample}", "clustering_full_ok.txt")
    params:
        cluster_dir=os.path.join(config['trycycler_cluster_path'], "{sample}"),
        threads=16,
        log=config['log'],
        queue="Priority"
    shell:
        """
        ARG=$(cat {input.depth})
        if (( $(echo "$ARG > 20" | bc -l) )); then
            rm -rf {params.cluster_dir}
            conda activate trycycler_env
            trycycler cluster --assemblies {input.assembly}/*.fasta --threads {params.threads} --reads {input.fastq} --out_dir {params.cluster_dir}
            conda deactivate
        fi
        """


