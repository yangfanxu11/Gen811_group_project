#!/bin/bash

#####################################################################
# Script: 01_leafcutter_prep.sh
# Purpose: Prepare data for Leafcutter analysis
# - Convert STAR SJ.out.tab to Leafcutter format
# - Create sample list and group files
#####################################################################

set -e
set -u

echo "======================================================================"
echo "STEP 1: LEAFCUTTER PREPARATION"
echo "======================================================================"

# Define paths
PROJECT_DIR="/home/users/yx1040/Gen811_group_project"
BAM_DIR="${PROJECT_DIR}/data/RNASEQ/alignment/star"
OUTPUT_DIR="${PROJECT_DIR}/data/leafcutter"
JUNC_DIR="${OUTPUT_DIR}/juncfiles"

# Create output directories
mkdir -p ${JUNC_DIR}

# Define samples (excluding KD3 and KD4)
WT_SAMPLES="WT1 WT2 WT3 WT4"
KD_SAMPLES="KD1 KD2"

echo ""
echo "[1.1] Converting STAR SJ.out.tab to junc format..."
echo "----------------------------------------------------------------------"

# Convert each STAR junction file
for sample in ${WT_SAMPLES} ${KD_SAMPLES}; do
    echo "Processing: ${sample}"
    
    SJ_FILE="${BAM_DIR}/${sample}_SJ.out.tab"
    JUNC_FILE="${JUNC_DIR}/${sample}.junc"
    
    if [ ! -f "${SJ_FILE}" ]; then
        echo "ERROR: SJ file not found: ${SJ_FILE}"
        exit 1
    fi
    
    # Convert STAR SJ.out.tab to simple BED6 format
    # STAR columns: chr start end strand motif annotated uniq_reads multi_reads max_overhang
    # We need: chr start end name score strand
    # Only keep junctions with strand info (column 4: 1=+, 2=-)
    awk 'BEGIN {OFS="\t"} 
         $4 == 1 {print $1, $2-1, $3, "JUNC", $7, "+"}
         $4 == 2 {print $1, $2-1, $3, "JUNC", $7, "-"}' \
         ${SJ_FILE} > ${JUNC_FILE}
    
    echo "  ✓ ${JUNC_FILE}"
    echo "  Junctions: $(wc -l < ${JUNC_FILE})"
done

echo ""
echo "[1.2] Creating sample list..."
echo "----------------------------------------------------------------------"

JUNCFILE_LIST="${JUNC_DIR}/juncfile_list.txt"
> ${JUNCFILE_LIST}

for sample in ${WT_SAMPLES} ${KD_SAMPLES}; do
    echo "${JUNC_DIR}/${sample}.junc" >> ${JUNCFILE_LIST}
done

echo "  ✓ ${JUNCFILE_LIST}"
echo "  Samples: $(wc -l < ${JUNCFILE_LIST})"

echo ""
echo "[1.3] Creating groups file..."
echo "----------------------------------------------------------------------"

GROUPS_FILE="${JUNC_DIR}/groups_file.txt"
> ${GROUPS_FILE}

for sample in ${WT_SAMPLES}; do
    echo "${sample} WT" >> ${GROUPS_FILE}
done

for sample in ${KD_SAMPLES}; do
    echo "${sample} KD" >> ${GROUPS_FILE}
done

echo "  ✓ ${GROUPS_FILE}"
cat ${GROUPS_FILE}

echo ""
echo "======================================================================"
echo "STEP 1 COMPLETE"
echo "======================================================================"
echo "Output: ${JUNC_DIR}"
echo "Ready for Leafcutter clustering"
echo "======================================================================"
