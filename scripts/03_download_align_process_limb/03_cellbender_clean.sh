#get Seurat.sif file to allow for creation of h5ad from count matrices
apptainer build Seurat.sif docker://satijalab/seurat:v4.3.0
#make h5ad for Cellbender
apptainer  exec  ${ScriptDIR}Seurat.sif Rscript Rscript.contra.intact.save.R

#now run cellbender
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
