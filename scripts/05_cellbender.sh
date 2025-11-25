#!/bin/bash

#SBATCH --job-name=cellbender         # Job name
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --time=                # Time limit (7 days)
#SBATCH -o %x_%j.out                  # Standard output (%x=job-name, %j=job-id)
#SBATCH -e %x_%j.err                  # Standard error (%x=job-name, %j=job-id)
#SBATCH --mail-type=ALL               # Email notifications
#SBATCH --mail-user=  # Email for notifications
#SBATCH -A                 # Project allocation
#SBATCH -p                    # Partition to use (GPU node)
#SBATCH --ntasks-per-node=16
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=5300            # Memory per CPU

pwd; hostname; date

# Accept input parameters
INPUT_DIR=$1
OUTPUT_DIR=$2
FPR=$3
EPOCHS=$4

# Default values if not provided
: ${FPR:=0.01}
: ${EPOCHS:=300}

# Check for .h5ad and .h5 files in the input directory
mapfile -t H5_FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f \( -name "*.h5" -o -name "*.h5ad" \))

if [ ${#H5_FILES[@]} -gt 0 ]; then
    if [ ${#H5_FILES[@]} -gt 1 ]; then
        echo "Multiple h5/h5ad files found. Running in parallel..."
        for file in "${H5_FILES[@]}"; do
            BASENAME=$(basename "$file" .h5)   
            BASENAME=${BASENAME%.h5ad}         
            
            FILE_OUTPUT_DIR="$OUTPUT_DIR/$BASENAME"
            mkdir -p "$FILE_OUTPUT_DIR"
            
            sbatch --job-name=cellbender --nodes=1 --time=7-00:00 -o %x_%j.out -e %x_%j.err --mail-type=ALL --mail-user=ni2651va-s@student.lu.se -A lu2024-2-39 -p gpua100i --ntasks-per-node=16 --gres=gpu:1 --mem-per-cpu=5300 --wrap="cd '$FILE_OUTPUT_DIR' && cellbender remove-background --cuda --input '$file' --output '${BASENAME}.h5ad' --fpr '$FPR' --epochs '$EPOCHS'"
        done
    else
        echo "Single h5/h5ad file found. Running normally..."
        BASENAME=$(basename "${H5_FILES[0]}" .h5)
        BASENAME=${BASENAME%.h5ad}
        
        FILE_OUTPUT_DIR="$OUTPUT_DIR/$BASENAME"
        mkdir -p "$FILE_OUTPUT_DIR"
        cd "$FILE_OUTPUT_DIR" || exit 1

        cellbender remove-background --cuda --input "${H5_FILES[0]}" --output "${BASENAME}_filtered.h5ad" --fpr "$FPR" --epochs "$EPOCHS"
    fi
else
    echo "No h5 or h5ad files found in the input directory. Exiting."
    exit 1
fi
