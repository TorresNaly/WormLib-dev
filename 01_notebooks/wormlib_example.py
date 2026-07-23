#!/usr/bin/env python
# coding: utf-8
"""
WormLib Example Pipeline
========================
Thin script that imports and calls library functions from src/wormlib.py.
No duplicate function definitions — wormlib.py is the single source of truth.
"""

# ============================================================================
# 0. Imports
# ============================================================================
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tifffile as tiff

# Reporting and Utils
import csv
from datetime import datetime
from PIL import Image
from reportlab.pdfgen import canvas
from reportlab.platypus import Table, TableStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter

# Resolve src/ directory and add to path. Prefer the current working tree when
# launched from the repo root, then fall back to this script's location.
current_dir = Path().resolve()
script_dir = Path(__file__).resolve().parent
search_roots = [current_dir, *current_dir.parents, script_dir, *script_dir.parents]
src_dir = next((root / 'src' for root in search_roots if (root / 'src').is_dir()), None)
if src_dir is None:
    raise RuntimeError("Could not find the WormLib src directory.")
sys.path.insert(0, str(src_dir))
main_dir = Path(src_dir.parent)


from wormlib import (
    load_images, segmentation, get_cell_stage_and_size_filtered,
    classify_2cell, classify_4cell, embryo_segmentation, nuclear_segmentation,
    spot_detection, analyze_rna_density, line_scan,
    normalize_optional_channel,
)

try:
    from sklearn.exceptions import InconsistentVersionWarning
    warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
except ImportError:
    pass

print("Source directory:", src_dir)
print("Main directory:", main_dir)
print("Current working directory:", current_dir)

# ============================================================================
# 1. Configuration
# ============================================================================

# Paths
image_path = main_dir / "data/230713_Lp306_L4440_11/230713_Lp306_L4440_11_R3D.dv"
image_ref = main_dir / "data/230713_Lp306_L4440_11/230713_Lp306_L4440_11_R3D_REF.dv"
output_directory = current_dir / "output_temp"
output_directory.mkdir(parents=True, exist_ok=True)

# Model paths
models_dir = main_dir / "models"
model_2cell_path = models_dir / "2-cell_classification_RFmodel.joblib"
model_4cell_path = models_dir / "4-cell_classification_RFmodel.joblib"

# Microscope parameters
voxel_size = (1448, 450, 450)           # Z, Y, X in nm
spot_radius_ch0 = (1409, 340, 340)      # PSF for mCherry channel
spot_radius_ch1 = (1283, 310, 310)      # PSF for Cy5 channel

# Channel names (set to None if the channel does not exist)
channel_names = {
    'Cy5': "mRNA1",
    'mCherry': "mRNA2",
    'FITC': None,
    'DAPI': "DAPI",
    'brightfield': "brightfield",
}
channel_indices = {
    'Cy5': 0,
    'mCherry': 1,
    'FITC': 2,
    'DAPI': 3,
    'brightfield': None,
}
Cy5 = normalize_optional_channel(channel_names.get('Cy5'))
mCherry = normalize_optional_channel(channel_names.get('mCherry'))

# Pipeline flags
run_cell_segmentation = True
run_cell_classifier = True
run_spot_detection = True
run_mRNA_heatmaps = True
run_rna_density_analysis = True
run_line_scan_analysis = True

# Segmentation parameters
embryo_diameter = 500
nuclei_diameter = 70

# ============================================================================
# 2. Load images
# ============================================================================

bf_result = None
if image_ref is not None:
    bf_result = load_images(
        image_path=str(image_ref),
        output_directory=output_directory,
        channel_names={'brightfield': channel_names.get('brightfield')},
        slice_to_plot=12,
    )

color_result = load_images(
    image_path=str(image_path),
    output_directory=output_directory,
    channel_names=channel_names,
    channel_indices=channel_indices,
    slice_to_plot=12,
)

# Unpack results
bf = bf_result['bf'] if bf_result else None
image_name = bf_result['image_name'] if bf_result else (color_result['image_name'] if color_result else "unknown")

image_Cy5 = color_result['image_Cy5'] if color_result else None
image_mCherry = color_result['image_mCherry'] if color_result else None
image_nuclei = color_result['image_nuclei'] if color_result else None
Cy5_array = color_result['Cy5_array'] if color_result else None
mCherry_array = color_result['mCherry_array'] if color_result else None
grid_width = color_result.get('grid_width', 80) if color_result else 80
grid_height = color_result.get('grid_height', 80) if color_result else 80


def first_analysis_shape():
    for array in (Cy5_array, mCherry_array):
        if array is not None:
            return array.shape[-2:]
    for image in (image_Cy5, image_mCherry, image_nuclei, bf):
        if image is not None:
            return image.shape[-2:]
    return None

# ============================================================================
# 3. Segmentation
# ============================================================================

masks_cytosol, masks_nuclei = None, None
features_df = None
df_long = None
cell_stage = "no-nuclei"
fallback_to_embryo = False

if run_cell_segmentation and bf is not None and image_nuclei is not None:
    image_cytosol = bf
    second_image_cytosol = image_nuclei
    masks_cytosol, masks_nuclei, _, _ = segmentation(
        image_cytosol, image_nuclei, second_image_cytosol,
        output_directory=output_directory,
    )
    segmentation_filename = os.path.join(output_directory, f'cell_segmentation_{image_name}.png')
    plt.savefig(segmentation_filename)
    plt.close()

    # Determine cell stage
    cell_stage, nuclei_sizes, masks_filtered = get_cell_stage_and_size_filtered(
        masks_nuclei, voxel_size,
    )

    # Run classifier
    if run_cell_classifier:
        classifier_attempted = False
        if cell_stage == "2-cell":
            classifier_attempted = True
            features_df = classify_2cell(
                masks_cytosol=masks_cytosol, bf=bf,
                image_name=image_name, output_directory=output_directory,
                model_path=str(model_2cell_path), verbose=True,
            )
        elif cell_stage == "4-cell":
            classifier_attempted = True
            features_df = classify_4cell(
                masks_cytosol=masks_cytosol, bf=bf,
                image_name=image_name, output_directory=output_directory,
                model_path=str(model_4cell_path), verbose=True,
            )
        else:
            print(f"Skipping classifier: Stage '{cell_stage}' is not supported.")
            fallback_to_embryo = True

        if classifier_attempted and features_df is None:
            fallback_to_embryo = True

    # If classifier returned None (count mismatch or unsupported stage), disable downstream
    if features_df is None:
        run_cell_classifier = False

elif run_cell_segmentation and bf is None and image_nuclei is not None:
    print("Running nuclear-only segmentation because no brightfield/reference image was loaded.")
    masks_nuclei = nuclear_segmentation(image_nuclei)
    masks_cytosol = masks_nuclei
    cell_stage = "nuclei-only"
    run_cell_classifier = False
    tiff.imwrite(os.path.join(output_directory, "masks_nuclei.tif"), masks_nuclei.astype(masks_nuclei.dtype))
    tiff.imwrite(os.path.join(output_directory, "masks_cytosol.tif"), masks_cytosol.astype(masks_cytosol.dtype))

# Fallback: whole-embryo segmentation
if (masks_cytosol is None or fallback_to_embryo) and bf is not None and image_nuclei is not None:
    print("Running fallback whole-embryo segmentation...")
    masks_cytosol, masks_nuclei, _, _ = embryo_segmentation(
        bf, image_nuclei, image_name, output_directory,
        embryo_diameter=embryo_diameter, nuclei_diameter=nuclei_diameter,
    )
    run_cell_classifier = False
    features_df = None

if run_spot_detection and masks_cytosol is None:
    analysis_shape = first_analysis_shape()
    if analysis_shape is not None:
        print("No segmentation mask available; using a whole-image mask for spot detection.")
        masks_cytosol = np.ones(analysis_shape, dtype=np.uint16)
        run_cell_classifier = False
        features_df = None

# ============================================================================
# 4. Spot detection
# ============================================================================

list_spots_in_each_cell_ch0, list_spots_in_each_cell_ch1 = [], []
spots_post_clustering_ch0, spots_post_clustering_ch1 = None, None

if run_spot_detection and masks_cytosol is not None:
    if Cy5 is not None and Cy5_array is not None:
        rna_channel = Cy5
        detection_color = "red"
        spots_post_clustering_ch0, clusters_ch0, list_spots_in_each_cell_ch0, _ = spot_detection(
            Cy5_array, voxel_size, spot_radius_ch0, masks_cytosol,
            image_name=image_name, rna_channel=rna_channel,
            detection_color=detection_color, output_directory=str(output_directory),
        )

    if mCherry is not None and mCherry_array is not None:
        rna_channel = mCherry
        detection_color = "blue"
        spots_post_clustering_ch1, clusters_ch1, list_spots_in_each_cell_ch1, _ = spot_detection(
            mCherry_array, voxel_size, spot_radius_ch1, masks_cytosol,
            image_name=image_name, rna_channel=rna_channel,
            detection_color=detection_color, output_directory=str(output_directory),
        )

# ============================================================================
# 4.2 Save mRNA counts
# ============================================================================

sum_spots_ch0 = sum(list_spots_in_each_cell_ch0) if list_spots_in_each_cell_ch0 else None
sum_spots_ch1 = sum(list_spots_in_each_cell_ch1) if list_spots_in_each_cell_ch1 else None

if any(x is not None for x in [sum_spots_ch0, sum_spots_ch1]):
    data_wide = {'Image ID': image_name}
    if Cy5 is not None and sum_spots_ch0 is not None:
        data_wide[f'{Cy5} total molecules'] = sum_spots_ch0
    if mCherry is not None and sum_spots_ch1 is not None:
        data_wide[f'{mCherry} total molecules'] = sum_spots_ch1
    df_quantification = pd.DataFrame([data_wide])
    quantification_output = os.path.join(output_directory, f'total_mRNA_counts_{image_name}.csv')
    df_quantification.to_csv(quantification_output, index=False)
    print("Saved wide CSV with total abundance:")
    print(df_quantification)

    num_regions = max(len(list_spots_in_each_cell_ch0), len(list_spots_in_each_cell_ch1))
    if num_regions > 0:
        rows_long = []
        for i in range(num_regions):
            row = {'Image ID': image_name, 'region_id': i + 1}
            if Cy5 is not None and list_spots_in_each_cell_ch0:
                row[Cy5] = list_spots_in_each_cell_ch0[i] if i < len(list_spots_in_each_cell_ch0) else 0
            if mCherry is not None and list_spots_in_each_cell_ch1:
                row[mCherry] = list_spots_in_each_cell_ch1[i] if i < len(list_spots_in_each_cell_ch1) else 0
            if features_df is not None:
                row['label'] = features_df.at[i, "highest_confidence_label"]
                row['confidence'] = round(features_df.at[i, "prediction_confidence"], 3)
            rows_long.append(row)
        df_long = pd.DataFrame(rows_long)
        output_prefix = "per_cell" if features_df is not None else "per_region"
        long_output = os.path.join(output_directory, f'{output_prefix}_mRNA_counts_{image_name}.csv')
        df_long.to_csv(long_output, index=False)
        print("Saved per-cell CSV:" if features_df is not None else "Saved per-region CSV:")
        print(df_long)

# ============================================================================
# 5. Spatial analysis of mRNA
# ============================================================================

# 5.1 Heatmaps
if run_mRNA_heatmaps and masks_cytosol is not None:
    def _create_heatmap(spots_x, spots_y, title_suffix, rna_max=None):
        img_width, img_height = masks_cytosol.shape[1], masks_cytosol.shape[0]
        cell_width = img_width / grid_width
        cell_height = img_height / grid_height
        grid = np.zeros((grid_height, grid_width), dtype=int)
        for x, y in zip(spots_x, spots_y):
            cx, cy = int(x / cell_width), int(y / cell_height)
            if 0 <= cx < grid_width and 0 <= cy < grid_height:
                grid[cy, cx] += 1
        if rna_max is not None:
            fig, axs = plt.subplots(1, 2, figsize=(8, 4))
            axs[0].imshow(rna_max, cmap='gray')
            axs[0].set_title(f"{title_suffix} max projection")
            axs[0].axis("off")
            im = axs[1].imshow(grid, cmap='CMRmap', interpolation='nearest')
            axs[1].set_title(f"{title_suffix} heatmap")
            axs[1].axis("off")
            cbar = fig.colorbar(im, ax=axs[1], fraction=0.046, pad=0.06)
            cbar.set_ticks([])
            plt.tight_layout()
            plt.savefig(os.path.join(output_directory, f"{title_suffix}_heatmap.png"), dpi=300, bbox_inches='tight')
            plt.show()
            plt.close()
        return grid

    if spots_post_clustering_ch0 is not None:
        _create_heatmap(spots_post_clustering_ch0[:, 2], spots_post_clustering_ch0[:, 1],
                        title_suffix=Cy5, rna_max=image_Cy5)
    if spots_post_clustering_ch1 is not None:
        _create_heatmap(spots_post_clustering_ch1[:, 2], spots_post_clustering_ch1[:, 1],
                        title_suffix=mCherry, rna_max=image_mCherry)

# 5.2 RNA density plots
if run_rna_density_analysis and masks_cytosol is not None:
    rna_names = [Cy5, mCherry]
    rna_images = [image_Cy5, image_mCherry]
    for mRNA_name, image in zip(rna_names, rna_images):
        if image is not None:
            analyze_rna_density(
                image=image, masks_cytosol=masks_cytosol, colormap='PiYG',
                mRNA_name=mRNA_name, image_name=image_name, output_directory=str(output_directory),
            )

# 5.3 Line scan
if run_line_scan_analysis and masks_cytosol is not None:
    rna_names = [Cy5, mCherry]
    rna_images = [image_Cy5, image_mCherry]
    for mRNA_name, image in zip(rna_names, rna_images):
        if image is not None:
            line_scan(
                image=image, masks_cytosol=masks_cytosol, colormap='PiYG',
                mRNA_name=mRNA_name, image_name=image_name, output_directory=str(output_directory),
                run_cell_classifier=run_cell_classifier, features_df=features_df, df_long=df_long,
            )

# ============================================================================
# 6. Export PDF report
# ============================================================================



output_pdf_path = os.path.join(output_directory, "report.pdf")

# Collect output files sorted by creation time
output_file_paths = []
for filename in os.listdir(output_directory):
    if filename.lower().endswith((".png", ".csv")):
        output_file_paths.append(os.path.join(output_directory, filename))
sorted_files = sorted(output_file_paths, key=lambda f: os.path.getctime(f))

run_time = datetime.now()
current_date = run_time.date()

margin = 32
c = canvas.Canvas(output_pdf_path, pagesize=letter)
c.setFont("Times-Roman", 16)
c.drawString(margin, 728, f"{image_name}")
c.setFont("Times-Roman", 14)
c.drawString(margin, 713, f"Report Generated: {current_date}")

origin = 700
padding = 20

for file in sorted_files:
    if file.endswith(".png"):
        with Image.open(file) as img:
            w, h = img.size
            aspect = w / h
        width = 145.6 * aspect
        if width > 548:
            height = 548 / aspect
            width = 548
        else:
            height = 145.6
        title = os.path.basename(file)

        if origin >= height + padding + 15:
            c.setFont("Times-Roman", 12)
            c.drawString(margin, origin - 15, title)
            c.drawImage(file, margin, origin - height - padding, width, height)
            origin -= height + padding + 10
        else:
            c.showPage()
            c.setFont("Times-Roman", 16)
            c.drawString(margin, 728, f"{image_name}")
            origin = 710
            c.setFont("Times-Roman", 12)
            c.drawString(margin, origin, title)
            c.drawImage(file, margin, origin - height - padding, width, height)
            origin -= height + padding + 20

    elif file.endswith(".csv"):
        with open(file, newline="") as csv_file:
            reader = csv.reader(csv_file)
            data = list(reader)
        table = Table(data)
        table.setStyle(TableStyle([
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("BACKGROUND", (0, 0), (-1, 0), colors.grey),
            ("ALIGN", (0, 0), (-1, -1), "CENTER"),
            ("FONTNAME", (0, 0), (-1, 0), "Times-Roman"),
            ("BOTTOMPADDING", (0, 0), (-1, 0), 12),
            ("BACKGROUND", (0, 1), (-1, -1), colors.beige),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.black),
        ]))
        table_width, table_height = table.wrapOn(c, 400, 600)
        origin -= table_height + padding
        table.drawOn(c, margin, origin)

c.save()
print(f"PDF report saved at {output_pdf_path}")
