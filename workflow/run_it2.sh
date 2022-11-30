#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH -o ${PWD}/snakemake.%j.out
#SBATCH -e ${PWD}/snakemake.%j.err


# conda activate snakemake
snakemake --profile profiles/biowulf --verbose -p --use-conda --jobs 400 --use-envmodules --latency-wait 120 -T 0 --configfile $1
