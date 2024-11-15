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

# Define the workflow
rule all:
    input:
        expand(os.path.join(config['kaiju_path'], "{sample}.png"),sample=config['ont_fastq_sample_names']),
        expand(os.path.join(config['kaiju_path'], "{sample}_taxids.tsv"),sample=config['ont_fastq_sample_names']),
        expand(os.path.join(config['kaiju_path'], "{sample}_taxonkit.tsv"),sample=config['ont_fastq_sample_names'])

#Rule to assign taxonomy to contigs using Kaiju
rule run_kaiju:
    input:
        fasta=os.path.join(config['flye_rough_path'], "{sample}","flye01", "assembly.fasta")
    output:
        kaiju_out=os.path.join(config['kaiju_path'], "{sample}.out"),
        taxids_tsv=os.path.join(config['kaiju_path'], "{sample}_taxids.tsv")
    params:
        kaiju_db_nodes=config['kaiju_db_nodes'],
        kaiju_db_fmi=config['kaiju_db']
    shell:
        """
        source activate kaiju_env
        kaiju -t {params.kaiju_db_nodes} -f {params.kaiju_db_fmi} -i {input.fasta} -o {output.kaiju_out}
        awk '{{print $3}}' {output.kaiju_out} > {output.taxids_tsv}
        conda deactivate
        """

# Rule to convert taxon IDs to full lineage using TaxonKit
rule run_taxonkit:
    input:
        kaiju_out=os.path.join(config['kaiju_path'], "{sample}.out")
    output:
        taxonkit=os.path.join(config['kaiju_path'], "{sample}_taxonkit.tsv")
    params:
        taxonkit_data_dir=config['taxonkit_db']
    shell:
        """
        source activate taxonkit_env
        cat {input.kaiju_out} | taxonkit lineage -c -i 3 --data-dir {params.taxonkit_data_dir} > {output.taxonkit}
        conda deactivate
        """

# Rule for plotting the results of Kaiju assignments
rule plot_contigs:
    input:
        taxonkit=os.path.join(config['kaiju_path'], "{sample}_taxonkit.tsv"),
        assembly=os.path.join(config['flye_rough_path'], "{sample}", "flye01", "assembly_info.txt")
    output:
        plot=os.path.join(config['kaiju_path'], "{sample}.png")
    params:
        colors="#432a74,#8c61fc,#00babe,#ff8c00,#c7c4fa,#bcebeb,#ffe6cc,#ff7f00,#984ea3,#ffff33"
    shell:
        """
        python scripts/plot_contigs.py {input.taxonkit} {input.assembly} {output.plot} "{params.colors}"
        """
