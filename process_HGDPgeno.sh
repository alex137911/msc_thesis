#!/bin/bash

#SBATCH --account=def-vmooser
#SBATCH --job-name=process_HGDPgeno
#SBATCH --output=process_HGDPgeno.out
#SBATCH --error=process_HGDPgeno.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=2G

# --------------------------------------------------------------
# Script efficiency (39708679)
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 2
# CPU Utilized: 1-11:52:45
# CPU Efficiency: 49.84% of 2-23:59:26 core-walltime
# Job Wall-clock time: 1-11:59:43
# Memory Utilized: 1.16 GB
# Memory Efficiency: 5.78% of 20.00 GB

# --------------------------------------------------------------
# Load necessary modules
module purge
module load StdEnv/2023
module load gcc/12.3
module load bcftools/1.19

# Human Genome Diversity Project with 1000 Genomes supplement
# Subset data for overlapping variants with QC CARTaGENE variants
IN_DIR="/lustre06/project/6061810/shared/HGDP_1KG/unfiltered_vcfs"
CAG_FILTERED_DIR="/lustre07/scratch/chanalex/CARTaGENE/QC"
OUT_DIR="/lustre07/scratch/chanalex/HGDP-1KG/QC/HGDP-CAG_Subset"
HGDP_ONLY_DIR="/lustre07/scratch/chanalex/HGDP-1KG/QC/HGDP-CAG_Subset"

mkdir -p "$OUT_DIR"

# Chromosomes to process
# CHROMOSOMES=($(seq -f "chr%g" 1 22))
CHROMOSOMES=("chr21")

# Remove CARTaGENE sample IDs (from HGDP-CAG subset when using bcftools isec)
HGDP_IDS="/lustre06/project/6050814/chanalex/msc_thesis/Data/HGDP/complete_HGDPsamples.txt"

# Loop through chromosomes
for CHR in "${CHROMOSOMES[@]}"; do
    # Input files
    HGDP_VCF="$IN_DIR/gnomad.genomes.v3.1.2.hgdp_tgp.${CHR}.vcf.bgz"
    CAG_FILTERED_VCF="$CAG_FILTERED_DIR/${CHR}_filtered.vcf.gz"
    OUTPUT_VCF="$OUT_DIR/${CHR}_HGDP_subset.vcf.gz"
    HGDP_ONLY_VCF="${HGDP_ONLY_DIR}/${CHR}_HGDP_only.vcf.gz"

    # Check that input files exist
    if [[ -f "$HGDP_VCF" && -f "$CAG_FILTERED_VCF" ]]; then
        echo "Processing $CHR..."

        # Index quality controlled CAG VCFs as needed
        if [[ ! -f "${CAG_FILTERED_VCF}.tbi" ]]; then
            echo "Indexing $CAG_FILTERED_VCF..."
            bcftools index -t "$CAG_FILTERED_VCF"
        fi

        # Use bcftools isec to compute intserections between VCFs
        # -n=2: only retain variants present in BOTH files
        bcftools isec -n=2 -p "$OUT_DIR" "$HGDP_VCF" "$CAG_FILTERED_VCF"

        # Compress the intersection files with bgzip
        # Necessary for bcftools merge
        echo "Compressing intersection files..."
        bgzip -f "$OUT_DIR/0000.vcf"
        bgzip -f "$OUT_DIR/0001.vcf"

        # Index the compressed intersection files
        # Also necessary for bcftools merge
        echo "Indexing compressed intersection files..."
        bcftools index -t "$OUT_DIR/0000.vcf.gz"
        bcftools index -t "$OUT_DIR/0001.vcf.gz"

        # Merge the intersection into a single VCF
        bcftools merge -Oz -o "$OUTPUT_VCF" "$OUT_DIR/0000.vcf.gz" "$OUT_DIR/0001.vcf.gz"

        # Index the output VCF
        bcftools index -t "$OUTPUT_VCF"

        # Subset the VCF to only include HGDP samples
        bcftools view \
            -S "$HGDP_IDS" \
            --force-samples \
            -Oz \
            -o "$HGDP_ONLY_VCF" \
            "$OUTPUT_VCF"

        # Index the new file
        bcftools index -t "$HGDP_ONLY_VCF"
        
        # Remove intermediate files
        rm -f "$OUT_DIR/0000.vcf.gz" "$OUT_DIR/0001.vcf.gz" "$OUT_DIR/0000.vcf.gz.tbi" "$OUT_DIR/0001.vcf.gz.tbi"

        echo "Subset for $CHR completed: $OUTPUT_VCF"
    
    else
        echo "Input files for $CHR not found. Skipping..."
    fi
done

echo "Processing complete."

# --------------------------------------------------------------
# Load necessary modules
module purge
module load StdEnv/2020
module load plink/1.9b_6.21-x86_64

# Quality control of HGDP subset 
# (i.e., variants present in both HGDP and CARTaGENE datasets)
IN_DIR="/lustre07/scratch/chanalex/HGDP-1KG/QC/HGDP-CAG_Subset"
OUT_DIR="/lustre07/scratch/chanalex/HGDP-1KG/QC"

# Remove samples which fail gnomAD QC hard filters
# Remove samples which were contaminated
# Two samples per line for PLINK
REMOVE_FILE="${IN_DIR}/HGDP_1KG.PassedQC-PLINK.txt"

# Set thresholds
MAX_SAMPLE_MISSINGNESS=0.05  # 5% genotype missingness per sample (remove samples with call rate < 95%)
MIN_MAF=0.01                 # Minor allele frequency > 1%
MAX_MISSINGNESS=0.05         # 5% genotype missingness per variant (remove variants missing in > 5% of samples)
HWE_PVAL=1e-6                # Variants which depart Hardy-Weinberg equilibrium (p-value < 1e-6)

# Loop through each chromosome
for CHR in "${CHROMOSOMES[@]}"; do
  INPUT_FILE="$IN_DIR/${CHR}_HGDP_subset.vcf.gz"
  OUTPUT_FILE="$OUT_DIR/${CHR}_HGDP_filtered.vcf.gz"
  LOG_FILE="$OUT_DIR/${CHR}_HGDP_processing.log"

  echo "Processing $INPUT_FILE" > "$LOG_FILE"

  # Step 1: Convert VCF to PLINK binary format
  plink --vcf "$HGDP_ONLY_VCF" \
        --remove "$REMOVE_FILE" \
        --const-fid \
        --make-bed \
        --out "${OUT_DIR}/${CHR}_temp" \
        --real-ref-alleles >> "$LOG_FILE" 2>&1
  
  # Step 2: Remove gnomAD failed samples and apply filters
  plink --bfile "${OUT_DIR}/${CHR}_temp" \
        --maf "$MIN_MAF" \
        --geno "$MAX_MISSINGNESS" \
        --mind "$MAX_SAMPLE_MISSINGNESS" \
        --hwe "$HWE_PVAL" midp \
        --recode vcf bgz \
        --out "${OUT_DIR}/${CHR}_HGDP_filtered" >> "$LOG_FILE" 2>&1

  # Validate the output
  if [[ ! -s "${OUT_DIR}/${CHR}_HGDP_filtered.vcf.gz" ]]; then
      echo "Error: Output file ${OUT_DIR}/${CHR}_HGDP_filtered.vcf.gz is empty or missing." >> "$LOG_FILE"
      continue
  fi

  # Consolidate PLINK logs into one file
  cat "${OUT_DIR}/${CHR}_temp.log" "${OUT_DIR}/${CHR}_HGDP_filtered.log" >> "$LOG_FILE"
  rm -f "${OUT_DIR}/${CHR}_temp.log" "${OUT_DIR}/${CHR}_HGDP_filtered.log" \
        "${OUT_DIR}/${CHR}_temp".*

  echo "Finished processing $INPUT_FILE" >> "$LOG_FILE"

  # Cleanup intermediate files
  rm -f "${OUT_DIR}/${CHR}_temp".* \
      "$HGDP_ONLY_VCF" \
      "${HGDP_ONLY_VCF}.tbi"
done