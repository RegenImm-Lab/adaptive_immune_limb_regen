import pandas as pd
import sys

def process_rowname(rowname):
    """Modify row names by splitting at '_' and taking the first part."""
    return rowname.split('_')[0]  

if len(sys.argv) != 3:
    print("Usage: python script.py input.tsv output.tsv")
    sys.exit(1)

# Get input and output file names
input_file = sys.argv[1]
output_file = sys.argv[2]

# Load the cell-by-gene matrix
df = pd.read_csv(input_file, sep='\t', index_col=0)

# Process row names
df.index = df.index.map(process_rowname)

# Save the processed matrix
df.to_csv(output_file, sep='\t')

print(f"Row names have been processed and saved to '{output_file}'.")
