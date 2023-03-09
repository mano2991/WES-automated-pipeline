#!/usr/bin/env bash

# Assign command-line arguments to variables
# $1: sample name, e.g. "sample_S1"
# $2: working directory
# $3: reference location
# $4: BED files
set -o errexit
set -o nounset

# Create output directory for variant calling results
mkdir -p output/"$1"/"$1"_Variant_calling

# Call variants using GATK HaplotypeCaller
docker run --rm \
  -v "$2":/data \
  -v "$3":/Reference \
  broadinstitute/gatk gatk HaplotypeCaller \
  --java-options "-Xmx100g" \
  -R /Reference/hg38.fasta \
  -I /data/output/"$1"/"$1"_Add_Read_Group/"$1"_alignment_RG.bam \
  -O /data/output/"$1"/"$1"_Variant_calling/"$1"_GATK.vcf.gz

# Extract variants of interest using BED files
docker run --rm \
  -v "$2":/data \
  -v "$3":/Reference \
  bioslimcontainers/tabix:1.7 tabix -h -R Reference/"$4" \
  /data/output/"$1"/"$1"_Variant_calling/"$1"_GATK.vcf.gz \
  > output/"$1"/"$1"_Variant_calling/"$1"_GATK_Covered.vcf

# Call variants using Mutect2
docker run --rm \
    -v "$2":/data \
    -v "$3":/Reference \
    broadinstitute/gatk Mutect2 \
    --java-options "-Xmx100g" \
    -R /Reference/hg38.fasta \
    -I /data/output/"$1"/"$1"_Add_Read_Group/"$1"_alignment_RG.bam \
    -O /data/output/"$1"/"$1"_Variant_calling/"$1"_GATK_Mutect2.vcf.gz
#Extract variants of interest using BED files
docker run --rm \
-v "$2":/data \
-v "$3":/Reference \
bioslimcontainers/tabix:1.7 \
tabix -h -R Reference/"$4" \
/data/output/"$1"/"$1"_Variant_calling/"$1"_GATK_Mutect2.vcf.gz \
> output/"$1"/"$1"_Variant_calling/"$1"_GATK_Mutect2_Covered.vcf
