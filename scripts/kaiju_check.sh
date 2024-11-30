#!/bin/sh
#SBATCH --error err_%x_%j
#SBATCH --output out_%x_%j
#SBATCH --job-name kaiju_check
#SBATCH --mail-user XX@YY.ZZ
#SBATCH --mail-type END,FAIL
#SBATCH -p Priority,Background
#SBATCH -N 1
#SBATCH --cpus-per-task=24

#assign taxonomy to contigs
mkdir -p kaiju_euk

for i in $(cat ont_samplenames.txt);do
    #assign taxonomy
    module load kaiju/1.10.1
    if [[ ! -f kaiju_euk/"$i".out ]]; then
	kaiju -t "/data/databases_food/Kaiju_nr_euk_20230510/nodes.dmp" -f "/data/databases_food/Kaiju_nr_euk_20230510/kaiju_db_nr_euk.fmi" -i 13_flyerough/"$i"/flye01/assembly.fasta -o kaiju_euk/"$i".out
	module unload kaiju/1.10.1
	#convert to text using taxonkit
	awk '{print $3}' kaiju_euk/"$i".out > kaiju_euk/"$i"_taxids.tsv
	module load taxonkit/0.16.0
	cat kaiju_euk/"$i".out | taxonkit lineage -c -i 3 --data-dir /data/databases/Taxonkit/ > kaiju_euk/"$i"_taxonkit.tsv
	module unload taxonkit/0.16.0
    fi
done
