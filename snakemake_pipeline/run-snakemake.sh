#!/bin/bash

mkdir -p log

snakemake \
   -s Snakefile \
   --quiet rules \
   --rerun-incomplete \
   --keep-going \
   --jobs 128 \
   --max-jobs-per-second 5 \
   --latency-wait 30 \
   --cluster-config cluster.yaml  \
   --cluster "sbatch \
              --partition=main \
              --output={cluster.log}_%j.out \
              --error={cluster.log}_%j.err \
              --account={username} \
              --job-name={cluster.name} \
              --time={cluster.time} \
              --cpus-per-task={cluster.cpus} \
              --mem={cluster.mem}"
