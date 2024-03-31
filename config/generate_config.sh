#!/bin/bash

#Generate config file for snakemake.

rm config.yaml
touch config.yaml

#1.fastq_sample_names:
cd fastq
ls *.fastq.gz | sed 's/^//; s/\.fastq.gz$//' > samplenames.txt
cp samplenames.txt ..
cd ..

#need samplenames.txt to exist first, then read into array
IFS=$'\n' read -d '' -r -a samples < samplenames.txt
# Join array elements with commas
joined_samples=$(printf ", %s" "${samples[@]}")
# Remove the leading comma and space
joined_samples="${joined_samples:2}"
# Write to the new file with desired format
echo "fastq_sample_names: [${joined_samples}]" >> config.yaml

#sequencing_summary.txt file
for file in "$PWD"/sequencing_summary*; do
    # Check if file exists and is a regular file
    if [ -f "$file" ]; then
        # Read the entire name of the file into a variable
        seq_summary="$file"
        break
    fi
done
echo "sequencing_summary: ${seq_summary}" >> config.yaml

#genome_sizes csv file
awk 'BEGIN { FS=","; printf "genome_sizes : {\n" } NR>1 { gsub(/^[[:space:]]*\"?/, "", $1); gsub(/[[:space:]]+$/, "", $3); printf "  \"%s\":\"%s\",\n", $1, $3 } END { printf "}\n" }' genome_sizes.csv >> config.yaml

#1. NanoPlot path
echo "nanoplot_path: ${PWD}/1_run_nanoplot" >> config.yaml
echo "nanoplot_threads: 8" >> config.yaml

#2. NanoComp path
echo "nanocomp_path: ${PWD}/2_run_nanocomp" >> config.yaml
echo "nanocomp_threads: 2" >> config.yaml

#3. Porechop_abi
echo "porechop_path: ${PWD}/3_porechop_abi" >> config.yaml
echo "porechop_threads: 32" >> config.yaml

#4. Filtlong
echo "filtlong_path: ${PWD}/4_filtlong" >> config.yaml
echo "filtlong_threads: 2" >> config.yaml

#5. Trycycler subsample
#issues with conda - installed with pip3 install git+https://github.com/rrwick/Trycycler.git
#follow with mamba install of dependencies (minimap2, miniasm)
echo "depth_path: ${PWD}/5_depth" >> config.yaml
echo "trycycler_subsample_path: ${PWD}/6_trycycler_subsample" >> config.yaml
echo "trycycler_subsample_threads: 12" >> config.yaml

#6. Assemblies
#Issues with flye install: mamba install -c bioconda -c conda-forge flye=2.9 python=3.9
echo "assemblies_work_path: ${PWD}/7_assemblies_work" >> config.yaml
echo "flye_threads: 12" >> config.yaml
echo "raven_threads: 8" >> config.yaml
echo "unicycler_threads: 32" >> config.yaml

#7. Clean Assemblies
echo "assemblies_path: ${PWD}/8_assemblies" >> config.yaml
echo "assemblers: [flye01, flye02, flye03, flyefull, raven01, raven02, raven03, ravenfull,unicycler01, unicycler02, unicycler03, unicyclerfull]" >> config.yaml

#8. Trcycler Clustering
echo "trycycler_cluster_path: ${PWD}/9_trycycler_cluster" >> config.yaml

#3. mash_path: 
#echo "mash_path: ${PWD}/3_run_mash" >> config.yaml

#5. mash_db
#echo "mash_db: /data/Food/analysis/R6564_NGS/Katie_F/databases/mash/refseq.genomes.k21s1000.msh" >> config.yaml

#6. spades_path
#echo "spades_path: ${PWD}/4_run_spades"  >> config.yaml

#7. logs
echo "log: ${PWD}/LOGS" >> config.yaml
