#!/bin/bash

#SBATCH --account=def-vmooser
#SBATCH --job-name=process_HGDPgeno
#SBATCH --output=process_HGDPgeno.out
#SBATCH --error=process_HGDPgeno.err
#SBATCH --time=2:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=2G

# --------------------------------------------------------------
# (39559386)


# --------------------------------------------------------------
# Load necessary modules
module purge
module load bcftools/1.19

# Human Genome Diversity Project with 1000 Genomes supplement
# Subset data for overlapping variants with QC CARTaGENE variants
IN_DIR="/lustre06/project/6061810/shared/HGDP_1KG/unfiltered_vcfs"
CAG_FILTERED_DIR="/lustre07/scratch/chanalex/CARTaGENE/QC"
OUT_DIR="/home/chanalex/scratch/HGDP-1KG/QC"

mkdir -p "$OUT_DIR"

# Chromosomes to process
# CHROMOSOMES=("chrX")
CHROMOSOMES=("chr22")

# Loop through chromosomes
for CHR in "${CHROMOSOMES[@]}"; do
    # Input files
    HGDP_VCF="$IN_DIR/gnomad.genomes.v3.1.2.hgdp_tgp.${CHR}.vcf.bgz.csi"
    CAG_FILTERED_VCF="$CAG_FILTERED_DIR/${CHR}_filtered.vcf.gz"
    OUTPUT_VCF="$OUT_DIR/${CHR}_HGDP_subset.vcf.gz"

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
        bcftools isec -@ 16 -n=2 -p "$OUT_DIR" "$HGDP_VCF" "$CAG_FILTERED_VCF"

        # Merge the intersection into a single VCF
        bcftools merge -Oz -o "$OUTPUT_VCF" "$OUT_DIR/0000.vcf" "$OUT_DIR/0001.vcf"

        # Index the output VCF
        bcftools index -t "$OUTPUT_VCF"
        
        # Remove intermediate files
        rm -f "$OUT_DIR/0000.vcf" "$OUT_DIR/0001.vcf"

        echo "Subset for $CHR completed: $OUTPUT_VCF"
    
    else
        echo "Input files for $CHR not found. Skipping..."
    fi
done

echo "Processing complete."

# --------------------------------------------------------------
# # Quality control of HGDP subset 
# # (i.e., variants present in both HGDP and CARTaGENE datasets)
# IN_DIR="/lustre06/project/6061810/shared/HGDP_1KG/unfiltered_vcfs"

# # Chromosomes to process
# # CHROMOSOMES=("chrX")
# CHROMOSOMES=("chr22")
# # CHROMOSOMES=($(seq -f "chr%g" 1 22) "chrX")

# # Set thresholds
# MAX_SAMPLE_MISSINGNESS=0.10  # 10% missingness per sample (genotype call rate > 90%)
# MIN_MAF=0.00                 # Minor allele frequency > 0%
# MAX_MISSINGNESS=0.10         # 10% missingness per variant (remove variants missing in > 10% of samples)
# HWE_PVAL=1e-15               # Hardy-Weinberg equilibrium p-value < 1e-15

# # Loop through each chromosome
# for CHR in "${CHROMOSOMES[@]}"; do
#   INPUT_FILE="$IN_DIR/gnomad.genomes.v3.1.2.hgdp_tgp.${CHR}.vcf.bgz"
#   OUTPUT_FILE="$OUT_DIR/${CHR}_filtered.vcf.bgz"
#   LOG_FILE="$OUT_DIR/${CHR}_processing.log"
#   {
#     echo "Processing $INPUT_FILE"

#     # Step 1: Filter for genotyping call rate > 90% (i.e., sample missingness per variant < 10%)
#     bcftools +fill-tags "$INPUT_FILE" -- -t F_MISSING | \
#       bcftools query -f '[%SAMPLE\t%F_MISSING\n]' | \
#       awk -v max_miss=${MAX_SAMPLE_MISSINGNESS} '$2 >= max_miss {print $1}' | \
#       bcftools view -S ^/dev/stdin "$INPUT_FILE" -Oz -o "${OUT_DIR}/${CHR}_GENO-CALL90.vcf.bgz"

#     # Log the number of samples removed
#     total_samples=$(bcftools query -l "$INPUT_FILE" | wc -l)
#     final_samples=$(bcftools query -l "${OUT_DIR}/${CHR}_GENO-CALL90.vcf.bgz" | wc -l)
#     echo "Step 1: Samples before filtering: $total_samples, after filtering (remove genotyping call rate < 90%): $final_samples, removed: $((total_samples - final_samples))" >> "$LOG_FILE"

#     # Step 2: Filter by minor allele frequency (MAF > 0%)
#     bcftools view -i "MAF > ${MIN_MAF}" "${OUT_DIR}/${CHR}_GENO-CALL90.vcf.bgz" -Oz -o "${OUT_DIR}/${CHR}_MAF_000.vcf.bgz"
#     before_variants=$(bcftools view -H "${OUT_DIR}/${CHR}_GENO-CALL90.vcf.bgz" | wc -l)
#     after_variants_step2=$(bcftools view -H "${OUT_DIR}/${CHR}_MAF_000.vcf.bgz" | wc -l)
#     echo "Step 2: Variants before filtering: $before_variants, after filtering (remove MAF < 1%): $after_variants_step2, removed: $((before_variants - after_variants_step2))" >> "$LOG_FILE"

#     # Step 3: Filter by missingness < 10% (variant missingness < 10%)
#     bcftools view -i "F_MISSING < ${MAX_MISSINGNESS}" "${OUT_DIR}/${CHR}_MAF_000.vcf.bgz" -Oz -o "${OUT_DIR}/${CHR}_MISS-10.vcf.bgz"
#     # bcftools view -i "F_MISSING < ${MAX_MISSINGNESS}" "${OUT_DIR}/${CHR}_GENO-CALL90.vcf.gz" -Oz -o "${OUT_DIR}/${CHR}_MISS-10.vcf.bgz"
#     after_variants_step3=$(bcftools view -H "${OUT_DIR}/${CHR}_MISS-10.vcf.bgz" | wc -l)
#     echo "Step 3: Variants before filtering: $after_variants_step2, after filtering (remove missingness > 10%): $after_variants_step3, removed: $((after_variants_step2 - after_variants_step3))" >> "$LOG_FILE"
    
#     # Step 4: Filter by Hardy-Weinberg Equilibrium (HWE p > 1e-15)
#     bcftools view -i "HWE > ${HWE_PVAL}" "${OUT_DIR}/${CHR}_MISS-10.vcf.bgz" -Oz -o "$OUTPUT_FILE"
#     after_variants_step4=$(bcftools view -H "$OUTPUT_FILE" | wc -l)
#     echo "Step 4: Variants before filtering: $after_variants_step3, after filtering (remove HWE < 1e-15): $after_variants_step4, removed: $((after_variants_step3 - after_variants_step4))" >> "$LOG_FILE"

#     # Validate final output
#     if [[ ! -s "$OUTPUT_FILE" ]]; then
#         echo "Error: Output file $OUTPUT_FILE is empty or missing." >> "$LOG_FILE"
#         continue
#     fi

#     # Log final sample count
#     final_samples=$(bcftools query -l "$OUTPUT_FILE" | wc -l)
#     echo "Final samples kept after filtering: $final_samples, removed: $((total_samples - final_samples))" >> "$LOG_FILE"

#     # Cleanup intermediate files
#     rm -f "${OUT_DIR}/${CHR}_GENO-CALL90.vcf.bgz" \
#           "${OUT_DIR}/${CHR}_MAF_000.vcf.bgz" \
#           "${OUT_DIR}/${CHR}_MISS-10.vcf.bgz" \
          
#     echo "Finished processing $INPUT_FILE" >> "$LOG_FILE"

#   } > "$LOG_FILE" 2>&1

# done