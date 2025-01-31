#!/usr/bin/bash

# Global variables
EXTRACTOR="H0"  # arbitrary name to save the results
EXTRACTOR_CLASS="tilingtool.extractors.dinov2_vit_giant.DinoV2ViTGiant"
EXTRACTOR_WEIGHTS="dinov2_bioptimus_g14"
EXTRACTOR_MEAN="[0.707223, 0.578729, 0.703617]"
EXTRACTOR_STD="[0.211883, 0.230117, 0.177517]"

#DATASET_ID="a3857b67-3f28-40e6-b1f9-4842a2ad49af"  # CHUV
DATASET_ID="3205eea0-5076-4294-8876-08cf98474603"  # UKER
ABSTRA_DATASET="/home/owkin/data/dataset-${DATASET_ID}"
HE_INPUT="${ABSTRA_DATASET}/HE/*eHnE.svs"
FOLDER_OUTPUT="/home/owkin/project/LUCAS_DATA/tiling-${DATASET_ID}/${EXTRACTOR}"
CONDA_ENV_TILINGTOOL="gbmerlangen"

eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV_TILINGTOOL

tilingtool extract \
        $HE_INPUT \
        $FOLDER_OUTPUT \
        --extractor $EXTRACTOR_CLASS \
        --extractor.init_args.weights $EXTRACTOR_WEIGHTS \
        --extractor.init_args.mixed_precision "true" \
        --extractor.init_args.compile "true" \
        --extractor.init_args.mean "$EXTRACTOR_MEAN" \
        --extractor.init_args.std "$EXTRACTOR_STD" \
        --matter-filter BRUNet \
        --num-workers 8
