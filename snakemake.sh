#!/bin/bash
#sbatch --cpus-per-task=2 --time=24:00:00 snakemake.sh

module load snakemake samtools macs bedtools sicer ucsc|| exit 1

sbcmd="sbatch --mem={cluster.mem} --cpus-per-task={threads} "
sbcmd+="--time={cluster.time}"
# --output={cluster.out}
snakemake -pr --keep-going --local-cores $SLURM_CPUS_PER_TASK \
    --jobs 50 --cluster-config cluster.json --cluster "$sbcmd" \
    --latency-wait 120 all
