#!/bin/sh
  for file in /home/anabwire/data/raw_data*.fastq
    do
    fastqc -f fastq --(no)extract -o /home/anabwire/results/fastqc_results ${file}
    done
