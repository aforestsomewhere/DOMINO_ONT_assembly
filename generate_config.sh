#!/bin/bash

#Generate config file for snakemake. Run from the main directory
rm config/config.yaml
touch config/config.yaml

#Default values for command line args
illumina_flag=false
ont_flag=false
shell_executable="/bin/bash"  # Default shell executable

# Function to show script usage
usage() {
    echo "Usage: $0 [--illumina] [--ont] [--help]"
    echo "  --illumina  Process Illumina data"
    echo "  --ont       Process Oxford Nanopore data"
    echo "  --shell_executable <path> Specify the shell executable (default: /bin/bash)"
    

    echo "  --help      Show this help message"
    exit 1
}
# Parse command-line options
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --illumina)
            illumina_flag=true
            shift # Move to the next argument
            ;;
        --ont)
            ont_flag=true
            shift # Move to the next argument
            ;;
        --shell_executable)
            # Ensure the next argument is provided and set the shell_executable path
            if [[ -n "$2" && ! "$2" =~ ^- ]]; then
                shell_executable="$2"
                shift 2  # Move past the flag and the argument
            else
                echo "Error: --shell_executable not provided."
                usage
            fi
            ;;
        --help)
            usage
            ;;
        *)
    echo "Unknown option: $1"
    usage
    ;;
    esac
done

# If neither flag is set, display a message and exit
if ! $illumina_flag && ! $ont_flag; then
    echo "No processing flag set. Use --illumina or --ont."
    usage
fi

#shell config
echo "shell_executable: $shell_executable" >> config/config.yaml
echo "shell_prefix: source /install/software/2024/minpy3/python3.10/bin/activate" >> config/config.yaml

# illumina_fastq_sample_names:
if $illumina_flag; then
    echo "Setting up Illumina data..."
    echo "raw_fastq_path: ${PWD}/illumina_fastq" >> config/config.yaml
    cd illumina_fastq
    ls *_R1.fastq.gz | sed 's/^//; s/\_R1.fastq.gz$//' > samplenames.txt
    cp samplenames.txt ..
    cd ..
    echo "Done."
fi

#ont_fastq_sample_names
if $ont_flag; then
    echo "Setting up ONT data..."
    echo "ont_fastq_path: ${PWD}/ont_fastq" >> config/config.yaml
    cd ont_fastq
    ls *.fastq.gz | sed 's/^//; s/\.fastq.gz$//' > ont_samplenames.txt
    cp ont_samplenames.txt ..
    cd ..
    echo "Done."
fi
#read genome sizes csv file
#genome_sizes csv file
#awk 'BEGIN { FS=","; printf "genome_sizes : {\n" } NR>1 { gsub(/^[[:space:]]*\"?/, "", $1); gsub(/[[:space:]]+$/, "", $3); printf "  \"%s\":\"%s\",\n", $1, $3 } END { printf "}\n" }' genome_sizes.csv >> config/config.yaml
#alt
awk 'BEGIN { FS=","; printf "genome_sizes : {\n" } NR>1 { gsub(/^[[:space:]]*"?/, "", $1); gsub(/[[:space:]]+$/, "", $3); printf "  \"%s\":\"%s\",\n", $1, $3 } END { printf "}\n" }' genome_sizes.csv >> config/config.yaml

if $illumina_flag; then
    #need samplenames.txt to exist first, then read into array
    IFS=$'\n' read -d '' -r -a samples < samplenames.txt
    # Join array elements with commas
    joined_samples=$(printf ", %s" "${samples[@]}")
    # Remove the leading comma and space
    joined_samples="${joined_samples:2}"
fi
#same for ont sample names
if $ont_flag; then
    IFS=$'\n' read -d '' -r -a ont_samples < ont_samplenames.txt
    # Join array elements with commas
    joined_ont_samples=$(printf ", %s" "${ont_samples[@]}")
    # Remove the leading comma and space
    joined_ont_samples="${joined_ont_samples:2}"
fi

# Write list of samplenames
echo "raw_fastq_sample_names: [${joined_samples}]" >> config/config.yaml
echo "ont_fastq_sample_names: [${joined_ont_samples}]" >> config/config.yaml
#2. fastqc-path
echo "fastqc_path: ${PWD}/1_run_fastqc" >> config/config.yaml
#3. multiqc path
echo "multiqc_path: ${PWD}/2_run_multiqc" >> config/config.yaml
#4. mash_path: 
echo "mash_path: ${PWD}/3_run_mash" >> config/config.yaml
#5. mash_db
echo "mash_db: /data/Food/analysis/R6564_NGS/Katie_F/databases/mash/refseq.genomes.k21s1000.msh" >> config/config.yaml
#6. jellyfish
echo "jellyfish_path: ${PWD}/4_run_jellyfish" >> config/config.yaml
# genomescope2
echo "genomescope2_path: ${PWD}/5_run_genomescope" >> config/config.yaml
#6. spades_path
echo "spades_path: ${PWD}/6_run_spades"  >> config/config.yaml
#7. bandage
echo "bandage_path: ${PWD}/7_run bandage" >> config/config.yaml
echo "bandage_binary: /data/Food/analysis/R6564_NGS/Katie_F/Bandage" >> config/config.yaml
#8. Prodigal
echo "faa_path: ${PWD}/8_run_prodigal"  >> config/config.yaml
#9. Porechop
echo "porechop_abi_path: ${PWD}/9_run_porechop_abi" >> config/config.yaml
echo "porechop_threads: 12" >> config/config.yaml
#10. Filtlong
echo "filtlong_path: ${PWD}/10_run_filtlong" >> config/config.yaml
echo "filtlong_threads: 2" >> config/config.yaml
#11. nanoplot
echo "nanoplot_path: ${PWD}/11_nanoplot" >> config/config.yaml
echo "nanoplot_threads: 8" >> config/config.yaml
#12. nanocomp
echo "nanocomp_path: ${PWD}/12_nanocomp" >> config/config.yaml
echo "nanocomp_threads: 2" >> config/config.yaml
#13. flye_rough
echo "flye_rough_path: ${PWD}/13_flyerough" >> config/config.yaml
echo "flye_rough_threads: 12" >> config/config.yaml
#logs
echo "log: ${PWD}/LOGS" >> config/config.yaml
#if depth filtering was already run, then add
