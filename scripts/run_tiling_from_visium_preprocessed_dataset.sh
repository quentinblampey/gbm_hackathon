#!/usr/bin/bash

# Adapted from https://github.com/owkin/visium/blob/lhe/updating_hipe_integration/workflow/rules/image_tiling.smk

# Global variables
EXTRACTOR="H0"  # arbitrary name to save the results
EXTRACTOR_CLASS="tilingtool.extractors.dinov2_vit_giant.DinoV2ViTGiant"
EXTRACTOR_WEIGHTS="dinov2_bioptimus_g14"
EXTRACTOR_MEAN="[0.707223, 0.578729, 0.703617]"
EXTRACTOR_STD="[0.211883, 0.230117, 0.177517]"

DATASET_ID="b5bd8840-cc6e-4de2-a78a-22315e9d1734" # CHUV
#DATASET_ID="40ce7e61-5ff6-467e-9972-0c85705f086c" # UKER
ABSTRA_DATASET="/home/owkin/data/dataset-${DATASET_ID}"
FOLDER_SPACE_RANGER_COUNT="${ABSTRA_DATASET}/spaceranger_count"
FOLDER_CONVERTED_HE="${ABSTRA_DATASET}/converted_he"
FOLDER_OUTPUT_BASE="/home/owkin/project/LUCAS_DATA/tiling-${DATASET_ID}"
FOLDER_OUTPUT_SPOT_TILING_CSV="${FOLDER_OUTPUT_BASE}/aligned_features"
FOLDER_OUTPUT_FEATURES="${FOLDER_OUTPUT_BASE}/aligned_features_${EXTRACTOR}"
CONDA_ENV_TILINGTOOL="gbmerlangen"

eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV_TILINGTOOL

mkdir -p "$FOLDER_OUTPUT_BASE"
mkdir -p "$FOLDER_OUTPUT_SPOT_TILING_CSV"
mkdir -p "$FOLDER_OUTPUT_FEATURES"
for sample_folder in "$FOLDER_SPACE_RANGER_COUNT"/*/; do
        if [ -d "$sample_folder" ]; then
                # Step 1: generate spot tiling csv
                sample_name=$(basename "$sample_folder");
                he_image="${FOLDER_CONVERTED_HE}/${sample_name}.tif"
                tissue_pos_csv="${FOLDER_SPACE_RANGER_COUNT}/${sample_name}/outs/spatial/tissue_positions.csv"
                spot_tiling_csv="${FOLDER_OUTPUT_SPOT_TILING_CSV}/${sample_name}_spot_tiling_reference.csv"
                # Grab spacerancer output folder from the tissue positions csv
                sr_out_folder=$(dirname $(dirname $tissue_pos_csv));

                touch "$spot_tiling_csv";
                echo "slide,space_ranger_out" > "$spot_tiling_csv";
                echo "$he_image,$sr_out_folder" >> "$spot_tiling_csv";


                # Step 2: run spot tiling
                features_path="${FOLDER_OUTPUT_FEATURES}/${sample_name}/features.npy"
                features_dir=$(dirname $features_path);

                tilingtool extract_visium \
                        $spot_tiling_csv \
                        $features_dir \
                        --extractor $EXTRACTOR_CLASS \
                        --extractor.init_args.weights $EXTRACTOR_WEIGHTS \
                        --extractor.init_args.mixed_precision "true" \
                        --extractor.init_args.compile "true" \
                        --extractor.init_args.mean "$EXTRACTOR_MEAN" \
                        --extractor.init_args.std "$EXTRACTOR_STD" \
                        --matter-filter BRUNet \
                        --num-workers 8 \
                        --use-tilingtool-matter-filter false

                # Copy the spot tiling csv to the features dir because it references the spaceranger folder
                # Might be helpful for traceability.
                cp $spot_tiling_csv $features_dir/spot_tiling_reference.csv

                # Move the content of the features folder to the parent folder
                mv $features_dir/${sample_name}/* $features_dir
                rm -rf $features_dir/${sample_name}
        fi
done
