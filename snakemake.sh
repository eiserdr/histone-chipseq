#!/bin/bash
#sbatch --cpus-per-task=1 --time=24:00:00 snakemake.sh

mkdir ../log

#specific versions:
module load snakemake/5.7.4 samtools/1.9 macs/2.2.5 bedtools/2.29.0 sicer/2-1.0.0 ucsc/392 idr/2.0.3 R|| exit 1
#module load snakemake samtools macs bedtools sicer ucsc idr R|| exit 1
###If you need to rerun snakemake, you probably have to run 'snakemake --unlock'
snakemake --unlock
###to run more cleanly, remove the comment to include --out, but you will lose information if there is an error.
sbcmd="sbatch --mem={cluster.mem} --cpus-per-task={threads} "
sbcmd+="--time={cluster.time} " # --out={cluster.out} "

snakemake -pr --keep-going --local-cores $SLURM_CPUS_PER_TASK \
    --jobs 200 --cluster-config cluster.json --cluster "$sbcmd" \
    --latency-wait 120 all
