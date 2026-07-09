#!/bin/bash
#SBATCH --account=csu95_alpine1
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=01:30:00
#SBATCH --partition=amilan
#SBATCH --output="${output_directory}/%j-report.out"

# Set the base directory and initialize paths
local_download_directory="/pl/active/onishimura_lab/PROJECTS/naly/bigfish/nd2-nikon_files/02_other-nd2_test_files"
input_directory="${local_download_directory}/input"
subdirectories=("$input_directory"/*)
python_script_path="/pl/active/onishimura_lab/PROJECTS/naly/bigfish/01_Ce-bigIFSH_single-embryo-6.py"

# Ensure SLURM_ARRAY_TASK_ID is within bounds of subdirectories array
export FOLDER_NAME="${subdirectories[${SLURM_ARRAY_TASK_ID}]}"
export OUTPUT_DIRECTORY="${local_download_directory}/output"

# Select image type
export DV_IMAGES="False"
export ND2_IMAGES="True"
export TIFF_IMAGES="False"

# Set channel names
export Cy5="set-3_mRNA"
export mCherry="lin-41_mRNA"
export FITC="nothing"
export DAPI="DAPI"
export brightfield="brightfield"

# Microscope parameters
export wavelength_cy5="670"
export wavelength_mCherry="610"
export na="1.42"
export refractive_index_medium="1.515"

# PSF parameters
export SPOT_RADIUS_CH0="1409,340,340"
export SPOT_RADIUS_CH1="1283,310,310"
export VOXEL_SIZE="1448,450,450"

# Feature selection
export PSF_CALCULATOR="True"
export SEGMENTATION="True"
export SPOT_DETECTION="True"

# Heatmaps
export RUN_mRNA_HEATMAPS="True"
export RUN_PROTEIN_HEATMAPS="True"
export ANALYZE_RNA_DENSITY="False"

# Mask generation
export GENERATE_DONUT_MASK="False"
export GENERATE_PGRANULE_MASK="False"

# Colocalization calculations
export CALCULATE_MEMBRANE_COLOCALIZATION="False"
export CALCULATE_NUCLEI_COLOCALIZATION="False"
export CALCULATE_PGRANULE_COLOCALIZATION="False"
export CALCULATE_mRNA_mRNA_COLOCALIZATION="False"

# Ensure SLURM_ARRAY_TASK_ID is within bounds of subdirectories array
folder_name="${subdirectories[${SLURM_ARRAY_TASK_ID}]}"

# Check if folder_name is valid
if [[ -d "$folder_name" ]]; then
    echo "Processing folder: $folder_name"
    echo "Files in folder:"
    ls -lh "$folder_name"  # List files with details for sanity check

    # Create a unique output directory for the current folder
    output_directory="${local_download_directory}/output/$(basename "$folder_name")"
    
    # Create output folder if it doesn't exist
    mkdir -p "$output_directory"

    # Set environment variables for the Python script execution
    export FOLDER_NAME="$folder_name"
    export OUTPUT_DIRECTORY="$output_directory"

    # Execute the Python script
    python "$python_script_path" "$folder_name" "$output_directory"

else
    echo "Invalid or missing directory: $folder_name"
fi

# Wait for all SLURM array tasks to complete before combining CSV files
if [[ ${SLURM_ARRAY_TASK_ID} -eq $((SLURM_ARRAY_TASK_COUNT - 1)) ]]; then
    wait
fi

# Combine all CSVs into a single file
combined_csv="${local_download_directory}/output/combined_quantification.csv"
experiment_csv=()

# Find all individual quantification CSVs in output subdirectories
for csv_file in "${local_download_directory}/output"/*/quantification_*.csv; do
    experiment_csv+=("$csv_file")
done

# Combine CSV files if data is found
if [[ ${#experiment_csv[@]} -gt 0 ]]; then
    # Use pandas in a short Python script to combine CSV files
    python - <<END
import pandas as pd
import glob
# Collect all quantification CSV files
csv_files = glob.glob("${local_download_directory}/output/*/quantification_*.csv")
combined_df = pd.concat([pd.read_csv(f) for f in csv_files], ignore_index=True)
combined_df.to_csv("${combined_csv}", index=False)
print("Combined CSV saved at ${combined_csv}")
END
else
    echo "No quantification data found to combine."
fi

