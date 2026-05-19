#!/bin/bash
set -e

echo "======================================================================"
echo "STEP 2: LEAFCUTTER DIFFERENTIAL SPLICING ANALYSIS"
echo "======================================================================"

# 使用绝对路径
PROJECT_DIR="/home/users/yx1040/Gen811_group_project"
JUNC_DIR="${PROJECT_DIR}/data/leafcutter/juncfiles"
OUT_DIR="${PROJECT_DIR}/data/leafcutter"
CLUSTER_SCRIPT="$HOME/leafcutter/clustering/leafcutter_cluster.py"
DIFF_SCRIPT="$HOME/leafcutter/scripts/leafcutter_ds.R"

echo ""
echo "[2.1] Clustering introns..."
echo "----------------------------------------------------------------------"
cd ${OUT_DIR}
python ${CLUSTER_SCRIPT} \
    -j ${JUNC_DIR}/juncfile_list.txt \
    -m 20 \
    -o leafcutter \
    -l 500000

echo "  ✓ Clustering complete"
echo "Total lines: $(zcat leafcutter_perind.counts.gz | wc -l)"

echo ""
echo "[2.2] Differential splicing analysis..."
echo "----------------------------------------------------------------------"
Rscript ${PROJECT_DIR}/leafcutter_ds_fixed.R \
    --num_threads 4 \
    --min_samples_per_intron 2 \
    --min_samples_per_group 2 \
    --min_coverage 10 \
    --output_prefix diff_splicing \
    leafcutter_perind.counts.gz \
    ${JUNC_DIR}/groups_file.txt

if [ -f "diff_splicing_cluster_significance.txt" ]; then
    cp diff_splicing_cluster_significance.txt diff_splicing_results.txt
    
    TOTAL=$(wc -l < diff_splicing_results.txt)
    echo "  ✓ Differential analysis complete"
    echo "  Total clusters tested: $((TOTAL - 1))"
else
    echo "  ✗ Error: diff_splicing_cluster_significance.txt not found"
    exit 1
fi

echo ""
echo "======================================================================"
echo "LEAFCUTTER ANALYSIS COMPLETE"
echo "======================================================================"
