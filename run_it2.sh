#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH -o ${PWD}/snakemake.%j.out
#SBATCH -e ${PWD}/snakemake.%j.err


module load singularity

mkdir -p TMP
export TMPDIR=TMP

# conda activate snakemake
snakemake --profile workflow/profiles/biowulf --verbose -p --use-conda --jobs 400 --default-resources "mem_mb=10000  " --use-envmodules --use-singularity --singularity-args " -B /vf,/spin1,/data,/fdb,/gpfs "  --latency-wait 120 -T 0 -s workflow/Snakefile --configfile $1
