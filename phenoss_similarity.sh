#!/bin/bash

# Check if both arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2

# Run the Python script
python ./src/similarity_score.py -input_dir "$INPUT_DIR" -output_dir "$OUTPUT_DIR"

# sbatch --cpus-per-task=1 --array=0-39 --mem=5G --time=5-00:00:00 --wrap="python ./src/similarity_score_local.py --input_file ... --output_dir ... --n_chunks 40 --index \$SLURM_ARRAY_TASK_ID"