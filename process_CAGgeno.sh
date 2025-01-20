#!/bin/bash

#SBATCH --account=def-vmooser
#SBATCH --job-name=process_CAGgeno
#SBATCH --output=process_CAGgeno.out
#SBATCH --error=process_CAGgeno.err
#SBATCH --time=6:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=2G

# --------------------------------------------------------------
# Script efficiency (37403891)
# State: COMPLETED (exit code 0)
# Nodes: 4
# Cores per node: 4
# CPU Utilized: 04:30:03
# CPU Efficiency: 6.87% of 2-17:32:16 core-walltime
# Job Wall-clock time: 04:05:46
# Memory Utilized: 28.67 MB
# Memory Efficiency: 0.09% of 32.00 GB

# --------------------------------------------------------------
# Load necessary modules
module purge
module load bcftools/1.19

# Input and output directories
IN_DIR="/lustre06/project/6061810/CERC_Private/Geno/CARTaGENE/Array/Pre-imputation"
OUT_DIR="/home/chanalex/scratch/CARTaGENE/QC"

mkdir -p "$OUT_DIR"

# Chromosomes to process
# CHROMOSOMES=("chrX")
CHROMOSOMES=($(seq -f "chr%g" 1 22) "chrX")

# Set thresholds
MAX_SAMPLE_MISSINGNESS=0.10  # 10% missingness per sample (genotype call rate > 90%)
MIN_MAF=0.00                 # Minor allele frequency > 0%
MAX_MISSINGNESS=0.10         # 10% missingness per variant (remove variants missing in > 10% of samples)
HWE_PVAL=1e-15               # Hardy-Weinberg equilibrium p-value < 1e-15

# Loop through each chromosome
for CHR in "${CHROMOSOMES[@]}"; do
  INPUT_FILE="$IN_DIR/${CHR}.CARTaGENEv1.1.array.vcf.gz"
  OUTPUT_FILE="$OUT_DIR/${CHR}_filtered.vcf.gz"
  LOG_FILE="$OUT_DIR/${CHR}_processing.log"
  {
    echo "Processing $INPUT_FILE"
    # Greer, P. J. et al. A reassessment of Hardy-Weinberg equilibrium filtering in large sample Genomic studies. 
    # 2024.02.07.24301951 Preprint at https://doi.org/10.1101/2024.02.07.24301951 (2024).

    # Step 1: Filter for genotyping call rate > 90% (i.e., sample missingness per variant < 10%)
    bcftools +fill-tags "$INPUT_FILE" -- -t F_MISSING | \
      bcftools query -f '[%SAMPLE\t%F_MISSING\n]' | \
      awk -v max_miss=${MAX_SAMPLE_MISSINGNESS} '$2 >= max_miss {print $1}' | \
      bcftools view -S ^/dev/stdin "$INPUT_FILE" -Oz -o "${OUT_DIR}/${CHR}_GENO-CALL90.vcf.gz"

    # Log the number of samples removed
    total_samples=$(bcftools query -l "$INPUT_FILE" | wc -l)
    final_samples=$(bcftools query -l "${OUT_DIR}/${CHR}_GENO-CALL90.vcf.gz" | wc -l)
    echo "Step 1: Samples before filtering: $total_samples, after filtering (remove genotyping call rate < 90%): $final_samples, removed: $((total_samples - final_samples))" >> "$LOG_FILE"

    # Step 2: Filter by minor allele frequency (MAF > 0%)
    bcftools view -i "MAF > ${MIN_MAF}" "${OUT_DIR}/${CHR}_GENO-CALL90.vcf.gz" -Oz -o "${OUT_DIR}/${CHR}_MAF_000.vcf.gz"
    before_variants=$(bcftools view -H "${OUT_DIR}/${CHR}_GENO-CALL90.vcf.gz" | wc -l)
    after_variants_step2=$(bcftools view -H "${OUT_DIR}/${CHR}_MAF_000.vcf.gz" | wc -l)
    echo "Step 2: Variants before filtering: $before_variants, after filtering (remove MAF < 1%): $after_variants_step2, removed: $((before_variants - after_variants_step2))" >> "$LOG_FILE"

    # Step 3: Filter by missingness < 10% (variant missingness < 10%)
    bcftools view -i "F_MISSING < ${MAX_MISSINGNESS}" "${OUT_DIR}/${CHR}_MAF_000.vcf.gz" -Oz -o "${OUT_DIR}/${CHR}_MISS-10.vcf.gz"
    after_variants_step3=$(bcftools view -H "${OUT_DIR}/${CHR}_MISS-10.vcf.gz" | wc -l)
    echo "Step 3: Variants before filtering: $after_variants_step2, after filtering (remove missingness > 10%): $after_variants_step3, removed: $((after_variants_step2 - after_variants_step3))" >> "$LOG_FILE"
    
    # Step 4: Filter by Hardy-Weinberg Equilibrium (HWE p > 1e-15)
    bcftools view -i "HWE > ${HWE_PVAL}" "${OUT_DIR}/${CHR}_MISS-10.vcf.gz" -Oz -o "$OUTPUT_FILE"
    after_variants_step4=$(bcftools view -H "$OUTPUT_FILE" | wc -l)
    echo "Step 4: Variants before filtering: $after_variants_step3, after filtering (remove HWE < 1e-15): $after_variants_step4, removed: $((after_variants_step3 - after_variants_step4))" >> "$LOG_FILE"

    # Validate final output
    if [[ ! -s "$OUTPUT_FILE" ]]; then
        echo "Error: Output file $OUTPUT_FILE is empty or missing." >> "$LOG_FILE"
        continue
    fi

    # Log final sample count
    final_samples=$(bcftools query -l "$OUTPUT_FILE" | wc -l)
    echo "Final samples kept after filtering: $final_samples, removed: $((total_samples - final_samples))" >> "$LOG_FILE"

    # Cleanup intermediate files
    rm -f "${OUT_DIR}/${CHR}_GENO-CALL90.vcf.gz" \
          "${OUT_DIR}/${CHR}_MAF_000.vcf.gz" \
          "${OUT_DIR}/${CHR}_MISS-10.vcf.gz" \
          
    echo "Finished processing $INPUT_FILE" >> "$LOG_FILE"

  } > "$LOG_FILE" 2>&1

done