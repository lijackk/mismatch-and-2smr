#!/bin/bash

mkdir -p log
snakemake \
   -s Snakefile_mvmr \
   --rerun-incomplete \
   --keep-going \
   --jobs 88 \
   --max-jobs-per-second 5 \
   --latency-wait 30 \
   --cluster-config cluster.yaml  \
   --cluster "sbatch \
              --output={cluster.log}_%j.out \
              --error={cluster.log}_%j.err \
              --job-name={cluster.name} \
              --time={cluster.time}  \
              --exclude=r6407,r6406  \
              --cpus-per-task={cluster.cpus}  \
              --mem={cluster.mem}"

