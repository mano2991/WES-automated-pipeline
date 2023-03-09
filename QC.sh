#!/usr/bin/env bash

# $1: sample name (e.g. "sample_S1")
# $2: work directory

# Set bash options to make the script safer and more predictable
set -o errexit  # Exit immediately if a command exits with a non-zero status
set -o nounset  # Treat unset variables as errors

# Create output directory and subdirectory for the current sample
mkdir -p output
mkdir -p output/${1}

# Run FastQC to check the quality of the input reads
mkdir -p output/${1}/QC_check
docker run --rm -v "$2":/data pegi3s/fastqc -t 30 /data/data/${1}_R1.fastq.gz -o /data/output/${1}/QC_check \
docker run --rm -v "$2":/data pegi3s/fastqc -t 30 /data/data/${1}_R2.fastq.gz -o /data/output/${1}/QC_check

# Use fastp to trim and filter the input reads
mkdir -p output/$1/trimedreads_fastp
docker run --rm -v "$2":/data biocontainers/fastp:v0.20.1_cv1 fastp -i data/${1}_R1.fastq.gz -I data/${1}_R2.fastq.gz \
  -o output/$1/trimedreads_fastp/${1}_1.fq.gz \
  -O output/$1/trimedreads_fastp/${1}_2.fq.gz \
  -j output/$1/trimedreads_fastp/${1}.json \
  -h output/$1/trimedreads_fastp/${1}.html \
  -w 16

# Run FastQC again to check the quality of the trimmed reads
mkdir -p output/${1}/Fastp_fastqc
docker run --rm -v"$2":/data pegi3s/fastqc -t 30 /data/output/$1/trimedreads_fastp/${1}_1.fq.gz -o /data/output/${1}/Fastp_fastqc \
docker run --rm -v"$2":/data pegi3s/fastqc -t 30 /data/output/$1/trimedreads_fastp/${1}_2.fq.gz -o /data/output/${1}/Fastp_fastqc
Note that the script now has spaces after the backslashes that indicate line continuation for improved readability.
