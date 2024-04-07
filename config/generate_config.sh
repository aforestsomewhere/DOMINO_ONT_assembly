#!/bin/bash

#Generate config file for snakemake.

rm config.yaml
touch config.yaml
mkdir -p LOGS

#shell config
echo "shell_executable: /bin/bash" >> config.yaml
echo "shell_prefix: source /install/software/anaconda3.6.b/bin/activate;" >> config.yaml
echo "raw_fastq_path: ${PWD}/fastq" >> config.yaml

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

#5. mash_path: 
echo "mash_path: ${PWD}/5_run_mash" >> config.yaml

#5a. mash_db
echo "mash_db: /data/Food/analysis/R6564_NGS/Katie_F/databases/mash/refseq.genomes.k21s1000.msh" >> config.yaml

#7. logs
echo "log: ${PWD}/LOGS" >> config.yaml
