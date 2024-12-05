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
        expand(os.path.join(config['kaiju_path'], "{sample}_taxonkit.tsv"),sample=config['ont_fastq_sample_names']),
	      os.path.join(config['kaiju_path'], "merge_kaiju_species.tsv"),
        expand(os.path.join(config['gtdb_path'],"in","{sample}.fasta"),sample=config['ont_fastq_sample_names']),
        os.path.join(config['gtdb_path'], "out","classify","gtdbtk.bac120.summary.tsv"),
	      os.path.join(config['quast_path'], "report.txt")


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
# Rule for generating merged table of kaiju results
rule kaiju_table:
    input:
        kaiju_out=expand(os.path.join(config['kaiju_path'], "{sample}.out"), sample=config['ont_fastq_sample_names'])
    output:
        kaiju_table=os.path.join(config['kaiju_path'], "merge_kaiju_species.tsv")
    params:
        kaiju_db_nodes=config['kaiju_db_nodes'],
        kaiju_db_names=config['kaiju_db_names']
    shell:
        """
        source activate kaiju_env
        kaiju2table -t {params.kaiju_db_nodes} -n {params.kaiju_db_names} -r species -o {output.kaiju_table} {input.kaiju_out}
        conda deactivate
        """
#Rule to prepare input folder for gtdb
rule prep_gtdb:
    input:
        assembly=os.path.join(config['flye_rough_path'], "{sample}","flye01", "assembly.fasta")
    output:
        assembly=os.path.join(config['gtdb_path'],"in","{sample}.fasta")
    params:
        queue = "Priority"
    shell:
        """
            mkdir -p $(dirname {output.assembly})
            cp {input.assembly} {output.assembly}
        """

#Rule to prepare input folder for gtdb
rule run_gtdb:
    input:
        assemblies=expand(os.path.join(config['gtdb_path'], "in", "{sample}.fasta"), sample=config['ont_fastq_sample_names'])
    output:
        summary = os.path.join(config['gtdb_path'], "out","classify","gtdbtk.bac120.summary.tsv"),
        all_output = directory(os.path.join(config['gtdb_path'],"out"))
    params:
        queue = "Priority",
        gtdb_db = config['GTDBTK_DATA_PATH'],
        outdir = os.path.join(config['gtdb_path'],"in")
    shell:
        """
            source activate gtdbtk-2.1.1
            export GTDBTK_DATA_PATH={params.gtdb_db}
            gtdbtk classify_wf --genome_dir $(dirname {input.assemblies[0]}) --out_dir {output.all_output} --skip_ani_screen -x fasta --cpus 24
            conda deactivate
        """
#Rule to run QUAST to generate assembly stats
rule run_quast:
    input:
        assemblies = expand(os.path.join(config['gtdb_path'], "in", "{sample}.fasta"), sample=config['ont_fastq_sample_names'])
    output:
        report = os.path.join(config['quast_path'], "report.txt")
    params:
        queue = "Priority",
	    outdir = os.path.join(config['quast_path'])
    shell:
        """
            source activate quast-5.3.0
            mkdir -p {params.outdir}
            quast -o {params.outdir} -L -t 8 {input.assemblies}
            conda deactivate
        """
