#!/usr/bin/env bash
#
# Run fixed-effect (FFX), random-effects (RFX), or both GWAMA models.
#
#   ./run_gwama.sh or meta_output cohort1.txt cohort2.txt [cohort3.txt ...]
#   ./run_gwama.sh or --model fixed  meta_output cohort1.txt cohort2.txt
#   ./run_gwama.sh or --model random meta_output cohort1.txt cohort2.txt
#   ./run_gwama.sh or --model both   meta_output cohort1.txt cohort2.txt
#
# Use qt instead of or for quantitative traits.
#
# Input files are produced by gwas_cohort_qc_with_gwama.R and already use
# GWAMA's default column names, so no --name_* flags are needed.

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
           With only a handful of cohorts, tau^2 is poorly estimated and the
           random-effects test is deflated: treat 'random' as a sensitivity
           analysis, not a primary result.

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

# Resolve GWAMA explicitly when the federated server keeps tools outside PATH.
# Resolution order: an exact FEDGX_GWAMA_BIN path, FEDGX_TOOLS_ROOT/GWAMA,
# then PATH.  serverSide_job.py sets FEDGX_TOOLS_ROOT for production jobs.
GWAMA_BIN="${FEDGX_GWAMA_BIN:-}"
if [[ -n "$GWAMA_BIN" ]]; then
    if [[ ! -x "$GWAMA_BIN" ]]; then
        echo "ERROR: FEDGX_GWAMA_BIN is not executable: $GWAMA_BIN" >&2
        exit 1
    fi
elif [[ -n "${FEDGX_TOOLS_ROOT:-}" && -x "${FEDGX_TOOLS_ROOT}/GWAMA" ]]; then
    GWAMA_BIN="${FEDGX_TOOLS_ROOT}/GWAMA"
elif command -v GWAMA >/dev/null 2>&1; then
    GWAMA_BIN="$(command -v GWAMA)"
else
    echo "ERROR: GWAMA was not found. Set FEDGX_GWAMA_BIN, place GWAMA at" >&2
    echo "       FEDGX_TOOLS_ROOT/GWAMA, or add it to PATH." >&2
    exit 1
fi

echo "Using GWAMA: $GWAMA_BIN"

# ---------------------------------------------------------------------------
# Validate the cohort files before handing anything to GWAMA
# ---------------------------------------------------------------------------
if [[ "$MODE" == "or" ]]; then
    REQUIRED_COLUMNS=(MARKERNAME EA NEA EAF N OR OR_95L OR_95U)
else
    REQUIRED_COLUMNS=(MARKERNAME EA NEA EAF N BETA SE)
fi

for file in "${COHORT_FILES[@]}"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: cohort file not found: $file" >&2
        exit 1
    fi

    if [[ "$file" == *.gz ]]; then
        echo "ERROR: GWAMA reads plain text; gunzip $file first." >&2
        exit 1
    fi

    header="$(head -1 "$file")"
    for column in "${REQUIRED_COLUMNS[@]}"; do
        if ! printf '%s\n' "$header" | tr '\t' '\n' | grep -qx -- "$column"; then
            echo "ERROR: $file is missing required column '$column' for mode '$MODE'." >&2
            echo "       header: $header" >&2
            exit 1
        fi
    done
done

OUTPUT_DIRECTORY="$(dirname -- "$OUTPUT_ROOT")"
mkdir -p -- "$OUTPUT_DIRECTORY"

# Sorted, so GWAMA's reference allele -- which it takes from the FIRST file in
# the list and applies to every SNP -- does not depend on the order the caller
# happened to pass the cohorts in. Without this, the sign of every beta can
# flip between otherwise identical runs.
FILELIST="${OUTPUT_ROOT}.gwama.in"
printf '%s\n' "${COHORT_FILES[@]}" | sort > "$FILELIST"

N_COHORTS="${#COHORT_FILES[@]}"

echo
echo "GWAMA input files (${N_COHORTS} cohorts, sorted):"
echo "-------------------------------------------"
cat "$FILELIST"
echo

if [[ "$MODEL" != "fixed" && "$N_COHORTS" -lt 5 ]]; then
    echo "NOTE: ${N_COHORTS} cohorts is few for a random-effects model;" >&2
    echo "      report fixed effects as the primary result." >&2
    echo >&2
fi

# ---------------------------------------------------------------------------
# Post-run checks
# ---------------------------------------------------------------------------
check_output() {
    # Reports how many variants were found in how many cohorts. If nothing
    # reaches n_studies == N_COHORTS, the marker names are not matching across
    # sites and no meta-analysis actually happened, even though GWAMA exits 0.
    local out_file="$1"
    local n_expected="$2"

    awk -v expected="$n_expected" -v f="$out_file" '
        NR == 1 {
            for (i = 1; i <= NF; i++) if ($i == "n_studies") col = i
            if (!col) { print "WARNING: no n_studies column in " f > "/dev/stderr"; exit 0 }
            next
        }
        col { counts[$col]++; if ($col + 0 > max) max = $col + 0 }
        END {
            if (!col) exit 0
            print "  n_studies distribution:"
            for (k = 1; k <= expected; k++)
                if (k in counts) printf "    in %d cohort(s): %d variants\n", k, counts[k]
            if (max < expected) {
                printf "ERROR: no variant present in all %d cohorts (max n_studies = %d).\n", expected, max > "/dev/stderr"
                printf "       Marker names are not matching across sites.\n" > "/dev/stderr"
                exit 1
            }
        }
    ' "$out_file"
}

report_log_warnings() {
    local log_file="$1"
    [[ -f "$log_file" ]] || return 0
    # GWAMA reports allele-frequency discrepancies above 30% to its log but
    # takes no action on them, so surface them here or they go unseen.
    local n
    n="$(grep -ci 'freq' "$log_file" || true)"
    if [[ "$n" -gt 0 ]]; then
        echo "  ${n} allele-frequency note(s) in ${log_file}:"
        grep -i 'freq' "$log_file" | head -10 | sed 's/^/    /'
        [[ "$n" -gt 10 ]] && echo "    ... ($((n - 10)) more)"
    fi
}

run_model() {
    local model="$1"
    local output_root="$2"
    local command=(
        "$GWAMA_BIN"
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

    check_output "${output_root}.out" "$N_COHORTS"
    report_log_warnings "${output_root}.log.out"
    echo
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

        COMPARE_SCRIPT="$SCRIPT_DIRECTORY/compare_gwama_models.py"
        if [[ ! -f "$COMPARE_SCRIPT" ]]; then
            echo "ERROR: --model both needs $COMPARE_SCRIPT, which is missing." >&2
            exit 1
        fi

        FIXED_ROOT="${OUTPUT_ROOT}.fixed"
        RANDOM_ROOT="${OUTPUT_ROOT}.random"
        COMPARISON_FILE="${OUTPUT_ROOT}.comparison.tsv"

        run_model fixed "$FIXED_ROOT"
        run_model random "$RANDOM_ROOT"

        python3 "$COMPARE_SCRIPT" \
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

echo "GWAMA completed successfully."
echo "Results:"
printf '  %s\n' "${RESULTS[@]}"
