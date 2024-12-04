#!/bin/bash

#Generate config file for snakemake. Run from the within the "isolates" folder

rm config/config.yaml
touch config/config.yaml

#Default values for command line args
illumina_flag=false
ont_flag=false
shell_executable="/bin/bash"  # Default shell executable
#shell_prefix="/install/software/2024/minpy3/python3.10/bin/activate" #default shell prefix

# Function to show usage
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
echo "config_file: ${PWD}/config/config.yaml" >> config/config.yaml
# illumina_fastq_sample_names:
if $illumina_flag; then
    echo "Setting up Illumina data..."
    echo "illumina_fastq_path: ${PWD}/illumina_fastq" >> config/config.yaml
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
echo "illumina_fastq_sample_names: [${joined_samples}]" >> config/config.yaml
echo "ont_fastq_sample_names: [${joined_ont_samples}]" >> config/config.yaml

#2. fastqc-path
echo "fastqc_path: ${PWD}/1_run_fastqc" >> config/config.yaml

#3. multiqc path
echo "multiqc_path: ${PWD}/2_run_multiqc" >> config/config.yaml

#4. mash_path: 
echo "mash_path: ${PWD}/3_run_mash" >> config/config.yaml

#5. mash_db
echo "mash_db: /data/databases_food/refseq.genomes.k21s1000.msh" >> config/config.yaml
# kaiju
echo "kaiju_path: ${PWD}/14_kaiju" >> config/config.yaml
echo "kaiju_db_nodes: /data/databases_food/Kaiju_nr_euk_20230510/nodes.dmp" >> config/config.yaml
echo "kaiju_db: /data/databases_food/Kaiju_nr_euk_20230510/kaiju_db_nr_euk.fmi" >> config/config.yaml
echo "kaiju_db_names: /data/databases_food/Kaiju_nr_euk_20230510/names.dmp" >> config/config.yaml
echo "taxonkit_db: /data/databases/Taxonkit" >> config/config.yaml
#6. jellyfish
echo "jellyfish_path: ${PWD}/4_run_jellyfish" >> config/config.yaml

# genomescope2
echo "genomescope2_path: ${PWD}/5_run_genomescope" >> config/config.yaml

#6. spades_path
echo "spades_path: ${PWD}/6_run_spades"  >> config/config.yaml
echo "spades_assemblies_path: ${PWD}/5_spades_assemblies" >> config/config.yaml
#7. bandage
echo "bandage_path: ${PWD}/7_run bandage" >> config/config.yaml
echo "bandage_binary: /data/Food/analysis/R6564_NGS/Katie_F/Bandage" >> config/config.yaml

#8. Prodigal
echo "faa_path: ${PWD}/8_run_prodigal"  >> config/config.yaml

#8. Porechop
echo "porechop_abi_path: ${PWD}/9_run_porechop_abi" >> config/config.yaml
echo "porechop_threads: 12" >> config/config.yaml

#9. Filtlong
echo "filtlong_path: ${PWD}/10_run_filtlong" >> config/config.yaml
echo "filtlong_threads: 2" >> config/config.yaml

#10. nanoplot
echo "nanoplot_path: ${PWD}/11_nanoplot" >> config/config.yaml
echo "nanoplot_threads: 8" >> config/config.yaml

#11. nanocomp
echo "nanocomp_path: ${PWD}/12_nanocomp" >> config/config.yaml
echo "nanocomp_threads: 2" >> config/config.yaml

#13. flye_rough
echo "flye_rough_path: ${PWD}/13_flyerough" >> config/config.yaml
echo "flye_rough_threads: 12" >> config/config.yaml

#12.subsampling
echo "subsamples_path: ${PWD}/12_subsamples" >> config/config.yaml

#13. assemblies
echo "assemblies_work_path: ${PWD}/13_assemblies" >> config/config.yaml
echo "assemblies_work_path_shallow: ${PWD}/13a_assemblies" >> config/config.yaml

#4. Clean Assemblies
echo "trycycler_path: ${PWD}/14_trycycler" >> config/config.yaml

#echo "assemblers: [flye01, flye02, flye03, flyefull, raven01, raven02, raven03, ravenfull,unicycler01, unicycler02, unicycler03, unicyclerfull]" >> config/config.yaml
echo "assemblers: [flye01, flye02, flye03, raven01, raven02, raven03, unicycler01, unicycler02, unicycler03]" >> config/config.yaml

#8. Trcycler Clustering
echo "trycycler_cluster_path: ${PWD}/15_trycycler_cluster" >> config/config.yaml

#18 Polishing
echo "polished_assemblies: ${PWD}/18_polishing" >> config/config.yaml
echo "pypolca_threads: 12" >> config/config.yaml
if [[! -f scripts/polypolish ]]; then
    #download latest binary release of polypolish and decompress
    curl -s https://api.github.com/repos/rrwick/Polypolish/releases/latest | grep "browser_download_url.*linux" | cut -d : -f 2,3 | tr -d \" | wget -qi -
    gunzip *polypolish-linux* && mv *polypolish-linux* scripts/
    tar -xvf scripts/*polypolish-linux* -C scripts/
    rm scripts/*.tar
    chmod u+x scripts/polypolish
    echo "polypolish_binary: ${PWD}/scripts/polypolish" >> config/config.yaml
fi

#8. gtdb_path
echo "gtdb_path: ${PWD}/15_gtdb" >> config/config.yaml
echo "GTDBTK_DATA_PATH: /data/databases_food/GTDB/release220/" >> config/config.yaml

echo "quast_path: ${PWD}/16_quast" >> config/config.yaml

#9. checkm_path
#echo "checkm2_path: ${PWD}/7_checkm2" >> config/config.yaml
#echo "checkm2_db: /data/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd" >> config/config.yaml
#echo "checkm2_binary: /data/Food/analysis/R6564_NGS/Katie_F/checkm2/bin/checkm2" >> config/config.yaml
#7. logs
echo "log: ${PWD}/LOGS" >> config/config.yaml
#if depth filtering was already run, then add
if [[ -f 13_depth/shallow.txt ]]; then
    ./workflow/assess_depth.sh
fi
