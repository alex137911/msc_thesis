#!/bin/bash

#SBATCH --account=def-vmooser
#SBATCH --job-name=ld_pruningCAG
#SBATCH --output=ld_pruningCAG.out
#SBATCH --error=ld_pruningCAG.err
#SBATCH --time=6:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=2G

# --------------------------------------------------------------
# Script efficiency (37555288)
# State: COMPLETED (exit code 0)
# Nodes: 4
# Cores per node: 4
# CPU Utilized: 00:08:53
# CPU Efficiency: 5.76% of 02:34:08 core-walltime
# Job Wall-clock time: 00:09:38
# Memory Utilized: 4.82 GB
# Memory Efficiency: 15.06% of 32.00 GB

# --------------------------------------------------------------
# Load necessary modules
module purge
module load StdEnv/2020
module load plink/1.9b_6.21-x86_64

# Directories
IN_DIR="/home/chanalex/scratch/CARTaGENE/QC"
OUT_DIR="/home/chanalex/scratch/CARTaGENE/LD_pruning"
EXCLUDE_FILE="/lustre06/project/6061810/shared/IBD_segments/grch38_genome_gap/genome_gap_hg38_and_MHC.bed"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Step 1: Convert VCFs to PLINK binary format
echo "Converting VCFs to PLINK format..."
for vcf in "$IN_DIR"/*_filtered.vcf.gz; do
    base=$(basename "$vcf" _filtered.vcf.gz)
    plink --vcf "$vcf" --make-bed --out "$OUT_DIR/$base"
done
echo "Conversion completed."

# Step 2: Combine all chromosomes into a single dataset
echo "Combining chromosomes into a single dataset..."
PLINK_FILES=$(ls "$OUT_DIR"/*.bed | sed 's/.bed//' > "$OUT_DIR/merge_list.txt")
plink --merge-list "$OUT_DIR/merge_list.txt" --make-bed --out "$OUT_DIR/merged_data"
echo "Chromosome merge completed."

# Step 3: Perform LD Pruning
echo "Performing LD pruning..."
plink --bfile "$OUT_DIR/merged_data" \
      --indep-pairwise 1000 50 0.05 \
      --exclude "$EXCLUDE_FILE" \
      --out "$OUT_DIR/pruned_data"
echo "LD pruning completed."

# Step 4: Create pruned dataset
echo "Creating pruned dataset..."
plink --bfile "$OUT_DIR/merged_data" \
      --extract "$OUT_DIR/pruned_data.prune.in" \
      --make-bed --out "$OUT_DIR/final_pruned_data"
echo "Pruned dataset created."

# Step 5: Clean up intermediate files
rm chr*.bed chr*.bim chr*.fam chr*.log chr*.nosex merge_list.txt
echo "LD pruning pipeline completed successfully."