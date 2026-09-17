#!/usr/bin/env bash

# ------------------------------------------------------------
# Run GWAMA meta-analysis
#
# Usage:
#
# Binary/case-control:
#   ./run_gwama.sh or meta_output \
#       cohort1.gwama.txt \
#       cohort2.gwama.txt \
#       cohort3.gwama.txt
#
# Quantitative:
#   ./run_gwama.sh qt meta_output \
#       cohort1.gwama.txt \
#       cohort2.gwama.txt \
#       cohort3.gwama.txt
# ------------------------------------------------------------

set -euo pipefail


# ------------------------------------------------------------
# Check arguments
# ------------------------------------------------------------

if [ "$#" -lt 3 ]; then
    echo "Usage:"
    echo "  $0 <or|qt> <output_root> <cohort1> <cohort2> [...]"
    exit 1
fi


MODE="$1"
OUTPUT_ROOT="$2"

shift 2

COHORT_FILES=("$@")


# ------------------------------------------------------------
# Check mode
# ------------------------------------------------------------

if [[ "$MODE" != "or" && "$MODE" != "qt" ]]; then
    echo "ERROR: mode must be 'or' or 'qt'"
    exit 1
fi


# ------------------------------------------------------------
# Check GWAMA
# ------------------------------------------------------------

if ! command -v GWAMA >/dev/null 2>&1; then
    echo "ERROR: GWAMA executable not found in PATH."
    exit 1
fi


# ------------------------------------------------------------
# Check cohort files
# ------------------------------------------------------------

for file in "${COHORT_FILES[@]}"; do

    if [[ ! -f "$file" ]]; then
        echo "ERROR: File not found: $file"
        exit 1
    fi

done


# ------------------------------------------------------------
# Create GWAMA input file list
# ------------------------------------------------------------

FILELIST="${OUTPUT_ROOT}.gwama.in"

> "$FILELIST"

for file in "${COHORT_FILES[@]}"; do
    echo "$file" >> "$FILELIST"
done


echo
echo "GWAMA file list:"
echo "----------------"
cat "$FILELIST"
echo


# ------------------------------------------------------------
# Run GWAMA
# ------------------------------------------------------------

if [[ "$MODE" == "or" ]]; then

    echo "Running binary/case-control GWAMA..."

    GWAMA \
        --filelist "$FILELIST" \
        --output "$OUTPUT_ROOT"

else

    echo "Running quantitative-trait GWAMA..."

    GWAMA \
        --filelist "$FILELIST" \
        --output "$OUTPUT_ROOT" \
        --quantitative

fi


# ------------------------------------------------------------
# Finished
# ------------------------------------------------------------

echo
echo "GWAMA completed."
echo
echo "Main result:"
echo "  ${OUTPUT_ROOT}.out"
echo
echo "Log files:"
echo "  gwama.log.out"
echo "  gwama.err.out"