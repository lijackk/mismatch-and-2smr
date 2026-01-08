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
              --account=lijack \
              --job-name={cluster.name} \
              --time={cluster.time} \
              --exclude=r6402,r6406,r6407,r6408,r6409,r6410 \
              --cpus-per-task={cluster.cpus} \
              --mem={cluster.mem}"