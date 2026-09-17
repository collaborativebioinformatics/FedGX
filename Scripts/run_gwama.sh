#!/usr/bin/env bash

# Run fixed-effect (FFX), random-effects (RFX), or both GWAMA models.
#
# Existing fixed-effect usage remains supported:
#   ./run_gwama.sh or meta_output cohort1.txt cohort2.txt [cohort3.txt ...]
#
# Select a model explicitly:
#   ./run_gwama.sh or --model fixed  meta_output cohort1.txt cohort2.txt
#   ./run_gwama.sh or --model random meta_output cohort1.txt cohort2.txt
#   ./run_gwama.sh or --model both   meta_output cohort1.txt cohort2.txt
#
# Use qt instead of or for quantitative traits.

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  run_gwama.sh <or|qt> [--model <fixed|random|both>] \
      <output_root> <cohort1> <cohort2> [cohort3 ...]

Arguments:
  or       Binary/case-control input (OR and confidence intervals).
  qt       Quantitative-trait input (BETA and SE).

Options:
  --model  Meta-analysis model. Default: fixed.

Outputs when --model both is used:
  <output_root>.fixed.out
  <output_root>.random.out
  <output_root>.comparison.tsv
EOF
}

if [[ $# -lt 4 ]]; then
    usage >&2
    exit 1
fi

MODE="$1"
shift

MODEL="fixed"
if [[ "${1:-}" == "--model" ]]; then
    if [[ $# -lt 2 ]]; then
        echo "ERROR: --model requires fixed, random, or both." >&2
        usage >&2
        exit 1
    fi
    MODEL="$2"
    shift 2
elif [[ "${1:-}" == --model=* ]]; then
    MODEL="${1#--model=}"
    shift
fi

if [[ $# -lt 3 ]]; then
    echo "ERROR: provide an output root and at least two cohort files." >&2
    usage >&2
    exit 1
fi

OUTPUT_ROOT="$1"
shift
COHORT_FILES=("$@")

if [[ "$MODE" != "or" && "$MODE" != "qt" ]]; then
    echo "ERROR: trait mode must be 'or' or 'qt'." >&2
    exit 1
fi

if [[ "$MODEL" != "fixed" && "$MODEL" != "random" && "$MODEL" != "both" ]]; then
    echo "ERROR: model must be 'fixed', 'random', or 'both'." >&2
    exit 1
fi

if ! command -v GWAMA >/dev/null 2>&1; then
    echo "ERROR: GWAMA executable not found in PATH." >&2
    exit 1
fi

for file in "${COHORT_FILES[@]}"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: cohort file not found: $file" >&2
        exit 1
    fi
done

OUTPUT_DIRECTORY="$(dirname -- "$OUTPUT_ROOT")"
mkdir -p -- "$OUTPUT_DIRECTORY"

FILELIST="${OUTPUT_ROOT}.gwama.in"
printf '%s\n' "${COHORT_FILES[@]}" > "$FILELIST"

echo
echo "GWAMA input files (${#COHORT_FILES[@]} cohorts):"
echo "-------------------------------------------"
cat "$FILELIST"
echo

run_model() {
    local model="$1"
    local output_root="$2"
    local command=(
        GWAMA
        --filelist "$FILELIST"
        --output "$output_root"
    )

    if [[ "$MODE" == "qt" ]]; then
        command+=(--quantitative)
    fi

    if [[ "$model" == "random" ]]; then
        command+=(--random)
    fi

    echo "Running GWAMA: trait=$MODE model=$model output=$output_root"
    "${command[@]}"

    if [[ ! -s "${output_root}.out" ]]; then
        echo "ERROR: GWAMA did not create ${output_root}.out" >&2
        exit 1
    fi
}

SCRIPT_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

case "$MODEL" in
    fixed)
        run_model fixed "$OUTPUT_ROOT"
        RESULTS=("${OUTPUT_ROOT}.out")
        ;;
    random)
        run_model random "$OUTPUT_ROOT"
        RESULTS=("${OUTPUT_ROOT}.out")
        ;;
    both)
        if ! command -v python3 >/dev/null 2>&1; then
            echo "ERROR: python3 is required to compare fixed and random results." >&2
            exit 1
        fi

        FIXED_ROOT="${OUTPUT_ROOT}.fixed"
        RANDOM_ROOT="${OUTPUT_ROOT}.random"
        COMPARISON_FILE="${OUTPUT_ROOT}.comparison.tsv"

        run_model fixed "$FIXED_ROOT"
        run_model random "$RANDOM_ROOT"

        python3 "$SCRIPT_DIRECTORY/compare_gwama_models.py" \
            --fixed "${FIXED_ROOT}.out" \
            --random "${RANDOM_ROOT}.out" \
            --output "$COMPARISON_FILE"

        RESULTS=(
            "${FIXED_ROOT}.out"
            "${RANDOM_ROOT}.out"
            "$COMPARISON_FILE"
        )
        ;;
esac

echo
echo "GWAMA completed successfully."
echo "Results:"
printf '  %s\n' "${RESULTS[@]}"
