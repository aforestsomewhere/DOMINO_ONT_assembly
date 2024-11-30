#!/bin/sh
#SBATCH --error err_%x_%j
#SBATCH --output out_%x_%j
#SBATCH --job-name mash_dist
#SBATCH --mail-user XX@YY.ZZ
#SBATCH --mail-type END,FAIL
#SBATCH --cpus-per-task=10
#SBATCH -p Priority
#SBATCH -N 1

source activate mash_test

#make list of input genomes
ls *.fasta | sed 's/^//; s/\.fasta$//' > samplenames.txt

#i by j
for i in $(cat samplenames.txt); do
    for j in $(cat samplenames.txt); do
	mash dist "$i".fasta "$j".fasta -p 10 > "$i"_"$j".out
    done
done
