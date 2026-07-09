#!/bin/bash
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=01:30:00
#SBATCH --mem=32G
#SBATCH --partition=amilan
#SBATCH --output="logs/%j-report.out"  # Initially save in a default logs directory

# Set the base directory and initialize paths
local_download_directory="$1"
input_directory="${local_download_directory}/input"
subdirectories=("$input_directory"/*)
python_script_path="./WormLib.py"

# 1. Define PSF parameters and channel names (set to None if the channel does not exist)
export SPOT_RADIUS_CH0="1409,340,340"
export SPOT_RADIUS_CH1="1283,310,310"
export VOXEL_SIZE="1448,450,450"

export ch0="set-3_mRNA"
export ch1="erm-1_mRNA"
export ch2="AIR2-GFP"
export ch3="DAPI"
export brightfield="brightfield"

# 2. Select image type
export DV_IMAGES="True"
export ND2_IMAGES="False"
export TIFF_IMAGES="False"

# 3. Segmentation
export RUN_EMBRYO_SEGMENTATION="False"
export EMBRYO_DIAMETER="375"
export NUCLEI_DIAMETER="70"
export RUN_CELL_SEGMENTATION="True"
export CELL_DIAMETER="250"
export RUN_CELL_CLASSIFIER="True"

# 4. Spot detection
export SPOT_DETECTION="True"

# 5. Spatial analysis of mRNA
export RUN_mRNA_HEATMAPS="True"
export RUN_RNA_DENSITY_ANALYSIS="False"
export RUN_LINE_SCAN_ANALYSIS="True"

# 6. Colocalization analysis
export CALCULATE_GFP_MASK_COLOCALIZATION="True"
export CALCULATE_NUCLEI_COLOCALIZATION="True"
export GENERATE_DONUT_MASK="True"
export CALCULATE_MEMBRANE_COLOCALIZATION="True"
export RUN_CONCENTRIC_LAYERS_ANALYSIS="True"
export GENERATE_PGRANULE_MASK="True"
export CALCULATE_PGRANULE_COLOCALIZATION="True"
export CALCULATE_mRNA_mRNA_COLOCALIZATION="True"
export RUN_PROTEIN_HEATMAPS="False"


# Ensure SLURM_ARRAY_TASK_ID is within bounds of subdirectories array
folder_name="${subdirectories[${SLURM_ARRAY_TASK_ID}]}"

# Check if folder_name is valid
if [[ -d "$folder_name" ]]; then
    echo "Files in folder:"
    ls -lh "$folder_name"  # List files with details for sanity check

    # Create a unique output directory for the current folder
    output_directory="${local_download_directory}/output/$(basename "$folder_name")"
    
    # Create output folder if it doesn't exist
    mkdir -p "$output_directory"

    # === Save script snapshots ===
    timestamp=$(date +%Y%m%d)
    cp "$python_script_path" "${output_directory}/${timestamp}_WormLib.py"
    cp "$0" "${output_directory}/${timestamp}_run-WormLib.sh"


    # Overwrite any existing log file
    : > "${output_directory}/log.out"  # truncates log.out if it exists

    # Redirect stdout and stderr to log file (overwrite)
    exec > >(tee "${output_directory}/log.out") 2>&1


    # Set environment variables for the Python script execution
    export FOLDER_NAME="$folder_name"
    export OUTPUT_DIRECTORY="$output_directory"

    # Execute the Python script
    python "$python_script_path" "$folder_name" "$output_directory"

else
    echo "Invalid or missing directory: $folder_name"
fi


# === Combine quantification_cell_counts CSVs ===
combined_counts_csv="${local_download_directory}/output/combined_cell_counts.csv"

python - <<END
import pandas as pd
import glob
import os

csv_files = glob.glob("${local_download_directory}/output/*/quantification_cell_*.csv")

if csv_files:
    combined_df = pd.concat([pd.read_csv(f) for f in csv_files], ignore_index=True)
    combined_df.to_csv("${combined_counts_csv}", index=False)
    print(f"Combined cell counts CSV saved at {os.path.abspath('${combined_counts_csv}')}")
else:
    print("No quantification_cell_counts CSVs found to combine.")
END
# else
#     echo "No single cell counts CSVs found to combine."
# fi

# # === Combine total_data CSVs ===
# combined_total_counts="${local_download_directory}/output/combined_total_counts.csv"
# total_counts_files=()

# for csv_file in "${local_download_directory}/output"/*/total_*.csv; do
#     total_counts_files+=("$csv_file")
# done

# if [[ ${#total_counts_files[@]} -gt 0 ]]; then
#     python - <<END
# import pandas as pd
# import glob
# csv_files = glob.glob("${local_download_directory}/output/*/total_*.csv")
# combined_df = pd.concat([pd.read_csv(f) for f in csv_files], ignore_index=True)
# combined_df.to_csv("${combined_total_counts}", index=False)
# print("Combined total counts data CSV saved at ${combined_total_counts}")
# END
# else
#     echo "No total counts CSVs found to combine."
# fi


# # === Combine line_scan_density_data CSVs ===
# combined_line_scan_counts="${local_download_directory}/output/combined_line_scan_counts.csv"
# combined_line_scan_files=()

# for csv_file in "${local_download_directory}/output"/*/*_line_density_data.csv; do
#     combined_line_scan_files+=("$csv_file")
# done

# if [[ ${#combined_line_scan_files[@]} -gt 0 ]]; then
#     python - <<END
# import pandas as pd
# import glob
# csv_files = glob.glob("${local_download_directory}/output/*/*_line_density_data.csv")
# combined_df = pd.concat([pd.read_csv(f) for f in csv_files], ignore_index=True)
# combined_df.to_csv("${combined_line_scan_counts}", index=False)
# print("Combined line scan data CSV saved at ${combined_line_scan_counts}")
# END
# else
#     echo "No line scan CSVs found to combine."
# fi
