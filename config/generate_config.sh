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

#1. NanoPlot path
echo "nanoplot_path: ${PWD}/1_run_nanoplot" >> config.yaml
echo "nanoplot_threads: 8" >> config.yaml

#2. NanoComp path
echo "nanocomp_path: ${PWD}/2_run_nanocomp" >> config.yaml
echo "nanocomp_threads: 2" >> config.yaml
#3. mash_path: 
#echo "mash_path: ${PWD}/3_run_mash" >> config.yaml

#5. mash_db
#echo "mash_db: /data/Food/analysis/R6564_NGS/Katie_F/databases/mash/refseq.genomes.k21s1000.msh" >> config.yaml

#6. spades_path
#echo "spades_path: ${PWD}/4_run_spades"  >> config.yaml

#7. logs
echo "log: ${PWD}/LOGS" >> config.yaml
