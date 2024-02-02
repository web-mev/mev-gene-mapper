#!/bin/bash

INPUT_FILE=$1
ORGANISM=$2
INITIAL_GENE_IDS=$3
TARGET_GENE_IDS=$4

if [ $ORGANISM = "Human" ]
then
    ORG_MAPPING="/opt/resources/human_mappings.tsv"
elif [ $ORGANISM = "Mouse" ]
then
    ORG_MAPPING="/opt/resources/mouse_mappings.tsv"
else
    echo "Not a valid organism choice." >&2
    exit 1;
fi

Rscript /usr/local/bin/map.R \
    -f $INPUT_FILE \
    -m $ORG_MAPPING \
    -i $INITIAL_GENE_IDS \
    -t $TARGET_GENE_IDS