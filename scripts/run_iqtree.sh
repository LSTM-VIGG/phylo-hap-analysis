#!/bin/bash

# ==============================================================================
#  Run IQ-TREE on multiple FASTA alignments with a specified outgroup.
# ==============================================================================

# --- Configuration ---
# Directory containing your input FASTA files
FASTA_DIR="fastas"

# Directory to store all IQ-TREE output files
OUTPUT_DIR="iqtree_results"

# The name of the outgroup sample as it appears in your FASTA headers
OUTGROUP= "VBS00259-4651STDY7017186_1" #"SAMN01760624_0"

# Number of bootstrap replicates for branch support
BOOTSTRAPS=1000

# Number of CPU cores to use. 'AUTO' lets IQ-TREE decide the optimal number.
THREADS="AUTO"
#MEM="256GB"
# --- End of Configuration ---

# Exit immediately if a command exits with a non-zero status.
set -e

# Create the output directory if it doesn't already exist
echo "Creating output directory: $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

# Find all .fasta files in the specified directory and loop through them
for alignment_file in "$FASTA_DIR"/*.fasta; do

    # Check if the file exists to avoid errors with empty directories
    if [ -f "$alignment_file" ]; then
        echo "--------------------------------------------------"
        echo "Starting IQ-TREE for: $alignment_file"
        echo "--------------------------------------------------"

        # Get the base name of the file without the extension (e.g., "vgsc_focal")
        # This will be used to create a unique prefix for the output files.
        prefix=$(basename "$alignment_file" .fasta)

        # Construct and run the IQ-TREE command
        # Using 'iqtree2' which is the common command for version 2+.
        # If your command is 'iqtree', change it below.
        iqtree2 \
            -s "$alignment_file" \
            -m MFP \
            -bb "$BOOTSTRAPS" \
            -T "$THREADS"

        echo "--------------------------------------------------"
        echo "Finished processing $alignment_file."
        echo "Results are prefixed with '$prefix' in the '$OUTPUT_DIR' directory."
        echo "--------------------------------------------------"
    fi
done

echo "All IQ-TREE jobs have completed."
