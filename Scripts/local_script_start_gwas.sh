#!/usr/bin/env bash
#
# Site-local GWAS driver.
#
#   1. resolve configuration (JSON file or environment)
#   2. prepare REGENIE-compatible phenotype / covariate files from LDAK output
#   3. run REGENIE step 1 and step 2
#   4. standardise to GWAMA format with gwas_cohort_qc_with_gwama.R
#
# The final artefact is ${OUTDIR}/${POPULATIONID}.GWAMA.txt, which is what the
# NVFLARE client streams back to the server.

set -euo pipefail

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
CONFIG_FILE="${HOME}/fedgx_config/config.json"

# ---------------------------------------------------------------------------
# Configuration: JSON file first, environment overrides it
# ---------------------------------------------------------------------------
if [ -f "${CONFIG_FILE}" ]; then
  while IFS='=' read -r key value; do
    case "${key}" in
      PROGRAM)       CFG_PROGRAM="${value}" ;;
      POPULATIONID)  CFG_POPULATIONID="${value}" ;;
      TRAITTYPE)     CFG_TRAITTYPE="${value}" ;;
      BUILD)         CFG_BUILD="${value}" ;;
    esac
  done < <(python3 - "${CONFIG_FILE}" <<'PY'
import json, shlex, sys
with open(sys.argv[1], encoding='utf-8') as f:
    cfg = json.load(f)
mapping = {
    "program": "PROGRAM",
    "populationid": "POPULATIONID",
    "traittype": "TRAITTYPE",
    "build": "BUILD",
}
for key, env_name in mapping.items():
    value = cfg.get(key, "")
    if value is None:
        value = ""
    print(f"{env_name}={shlex.quote(str(value))}")
PY
)
fi

PROGRAM="${PROGRAM:-${CFG_PROGRAM:-}}"
POPULATIONID="${POPULATIONID:-${CFG_POPULATIONID:-}}"
TRAITTYPE="${TRAITTYPE:-${CFG_TRAITTYPE:-binary}}"
BUILD="${BUILD:-${CFG_BUILD:-GRCh38}}"

if [ -z "${PROGRAM}" ]; then
  echo "ERROR: PROGRAM is not set (config key 'program' or env PROGRAM)." >&2
  exit 1
fi
if [ -z "${POPULATIONID}" ]; then
  echo "ERROR: POPULATIONID is not set (config key 'populationid' or env POPULATIONID)." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# Paths
#
# client.py passes BFILE / PHENO_FILE / COVAR_FILE / OUTDIR explicitly, taken
# from the site's datasets.json. When run by hand without them, fall back to
# the <DATAPATH>/<POPULATIONID>/<POPULATIONID>_* convention.
# ---------------------------------------------------------------------------
DATAPATH_DEFAULT="/home/ubuntu/data"

if [ -z "${DATAPATH:-}" ]; then
  DATAPATH="${DATAPATH_DEFAULT}"
else
  echo "DATAPATH is set to: ${DATAPATH}"
fi
export DATAPATH

POPULATIONPATH="${DATAPATH}/${POPULATIONID}"

BFILE="${BFILE:-${POPULATIONPATH}/${POPULATIONID}_geno}"
COVAR_RAW="${COVAR_FILE:-${POPULATIONPATH}/${POPULATIONID}_geno.covar}"
PHENO_RAW="${PHENO_FILE:-${POPULATIONPATH}/${POPULATIONID}_pheno.pheno}"
OUTDIR="${OUTDIR:-${POPULATIONPATH}/gwas_out}"

WORKDIR="${OUTDIR}/work"
mkdir -p "${WORKDIR}"

# Binaries: client.py resolves these from tools.json and passes them in.
# Fall back to whatever is on PATH.
TOOL_BIN="${TOOL_BIN:-}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"

QC_SCRIPT="${SCRIPT_DIR}/gwas_cohort_qc_with_gwama.R"

echo "========================================"
echo "Site GWAS"
echo "========================================"
echo "Population : ${POPULATIONID}"
echo "Program    : ${PROGRAM}"
echo "Trait type : ${TRAITTYPE}"
echo "Build      : ${BUILD}"
echo "Genotypes  : ${BFILE}"
echo "Binary     : ${TOOL_BIN:-<from PATH>}"
echo "Phenotype  : ${PHENO_RAW}"
echo "Covariates : ${COVAR_RAW}"
echo "Output dir : ${OUTDIR}"
echo "========================================"

for f in "${BFILE}.bed" "${BFILE}.bim" "${BFILE}.fam" "${COVAR_RAW}" "${PHENO_RAW}"; do
  if [ ! -f "${f}" ]; then
    echo "ERROR: required input not found: ${f}" >&2
    exit 1
  fi
done

if [ ! -f "${QC_SCRIPT}" ]; then
  echo "ERROR: QC script not found: ${QC_SCRIPT}" >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# Prepare REGENIE inputs
#
# LDAK writes headerless FID IID <value> files. REGENIE requires a header line,
# and for a binary trait expects 0 = control, 1 = case. LDAK/PLINK conventions
# vary between 0/1 and 1/2, so detect rather than assume.
# ---------------------------------------------------------------------------
PHENO="${WORKDIR}/pheno.txt"
COVAR="${WORKDIR}/covar.txt"
PHENO_NAME="Y1"

if [ "${TRAITTYPE}" = "binary" ]; then
  # Detect the coding from the values actually present, then map once.
  VALUES="$(awk '{ v[$3] } END { n = asorti(v, k); for (i = 1; i <= n; i++) printf "%s ", k[i] }' "${PHENO_RAW}" 2>/dev/null \
            || awk '{ v[$3] } END { for (k in v) printf "%s ", k }' "${PHENO_RAW}")"
  echo "Phenotype values present: ${VALUES}"

  if printf '%s' "${VALUES}" | grep -qE '(^| )2( |$)'; then
    echo "Detected PLINK/LDAK 1/2 coding; remapping to 0/1 for REGENIE."
    MAP_CASE=2
  else
    echo "Detected 0/1 coding; passing through."
    MAP_CASE=1
  fi

  awk -v OFS='\t' -v name="${PHENO_NAME}" -v case_val="${MAP_CASE}" '
    BEGIN { print "FID", "IID", name }
    {
      if ($3 == case_val)                    out = 1
      else if ($3 == 0 || $3 == 1 || $3 == 2) out = 0
      else                                    out = "NA"
      print $1, $2, out
    }
  ' "${PHENO_RAW}" > "${PHENO}"
else
  awk -v OFS='\t' -v name="${PHENO_NAME}" '
    BEGIN { print "FID", "IID", name }
    { print $1, $2, $3 }
  ' "${PHENO_RAW}" > "${PHENO}"
fi

echo "Case/control counts in ${PHENO}:"
awk 'NR > 1 { c[$3]++ } END { for (k in c) printf "  %s: %d\n", k, c[k] }' "${PHENO}"

# Covariates: header is FID IID C1 C2 ... for however many columns exist.
awk -v OFS='\t' '
  NR == 1 {
    printf "FID\tIID"
    for (i = 3; i <= NF; i++) printf "\tC%d", i - 2
    printf "\n"
  }
  { $1 = $1; print }          # rebuild the row so the separator is a tab too
' "${COVAR_RAW}" > "${COVAR}"

NCOVAR=$(( $(head -1 "${COVAR}" | awk '{print NF}') - 2 ))
echo "Covariates: ${NCOVAR}"

# ---------------------------------------------------------------------------
# Run the GWAS
# ---------------------------------------------------------------------------
STEP1_OUT="${WORKDIR}/step1"
STEP2_OUT="${WORKDIR}/step2"

case "${PROGRAM}" in

  regenie)
    echo ""
    echo "Running REGENIE workflow"
    REGENIE="${TOOL_BIN:-$(command -v regenie || true)}"
    if [ -z "${REGENIE}" ] || [ ! -x "${REGENIE}" ]; then
      echo "ERROR: regenie binary not found (TOOL_BIN='${TOOL_BIN:-}')" >&2
      exit 1
    fi
    echo "Using regenie: ${REGENIE}"

    REGENIE_TRAIT_FLAGS=()
    if [ "${TRAITTYPE}" = "binary" ]; then
      REGENIE_TRAIT_FLAGS=(--bt)
    fi

    # REGENIE step 1 fits the whole-genome model and is meant to run on a
    # pruned subset of common variants, not the full set. Use a supplied list
    # if there is one, otherwise take an evenly spaced subset across the
    # genome, which keeps every chromosome represented.
    STEP1_SNPS="${STEP1_SNPS:-20000}"
    STEP1_EXTRACT="${STEP1_EXTRACT:-}"

    if [ -z "${STEP1_EXTRACT}" ]; then
      STEP1_EXTRACT="${WORKDIR}/step1_variants.snps"
      TOTAL_SNPS=$(wc -l < "${BFILE}.bim")
      if [ "${TOTAL_SNPS}" -le "${STEP1_SNPS}" ]; then
        cut -f2 "${BFILE}.bim" > "${STEP1_EXTRACT}"
      else
        awk -v want="${STEP1_SNPS}" -v total="${TOTAL_SNPS}" '
          { acc += want; if (acc >= total) { acc -= total; print $2 } }
        ' "${BFILE}.bim" > "${STEP1_EXTRACT}"
      fi
      echo "Step 1 variant subset: $(wc -l < "${STEP1_EXTRACT}") of ${TOTAL_SNPS}"
    else
      echo "Step 1 variant subset: ${STEP1_EXTRACT} (supplied)"
    fi

    echo ""
    echo "--- REGENIE step 1 ---"
    "${REGENIE}" \
      --step 1 \
      --bed "${BFILE}" \
      --extract "${STEP1_EXTRACT}" \
      --covarFile "${COVAR}" \
      --phenoFile "${PHENO}" \
      ${REGENIE_TRAIT_FLAGS[@]+"${REGENIE_TRAIT_FLAGS[@]}"} \
      --bsize 1000 \
      --lowmem \
      --lowmem-prefix "${WORKDIR}/tmp_rg" \
      --threads "${THREADS:-4}" \
      --out "${STEP1_OUT}"

    if [ ! -f "${STEP1_OUT}_pred.list" ]; then
      echo "ERROR: step 1 did not produce ${STEP1_OUT}_pred.list" >&2
      exit 1
    fi

    echo ""
    echo "--- REGENIE step 2 ---"
    STEP2_FLAGS=()
    if [ "${TRAITTYPE}" = "binary" ]; then
      STEP2_FLAGS=(--bt --firth --approx --pThresh 0.05)
    fi

    "${REGENIE}" \
      --step 2 \
      --bed "${BFILE}" \
      --covarFile "${COVAR}" \
      --phenoFile "${PHENO}" \
      --pred "${STEP1_OUT}_pred.list" \
      ${STEP2_FLAGS[@]+"${STEP2_FLAGS[@]}"} \
      --bsize 400 \
      --threads "${THREADS:-4}" \
      --out "${STEP2_OUT}"

    SUMSTATS="${STEP2_OUT}_${PHENO_NAME}.regenie"
    COL_ARGS=(
      --col-chr CHROM
      --col-pos GENPOS
      --col-id ID
      --col-ea ALLELE1
      --col-nea ALLELE0
      --col-eaf A1FREQ
      --col-beta BETA
      --col-se SE
      --col-n N
      --col-log10p LOG10P
      --col-test TEST
      --filter-test TRUE
      --test-value ADD
    )
    ;;

  plink|gcta|saige)
    echo "PROGRAM value '${PROGRAM}' is not yet supported." >&2
    exit 1
    ;;

  *)
    echo "Unknown PROGRAM value: ${PROGRAM}" >&2
    exit 1
    ;;
esac

if [ ! -s "${SUMSTATS}" ]; then
  echo "ERROR: expected summary statistics not found: ${SUMSTATS}" >&2
  exit 1
fi

echo ""
echo "Summary statistics: ${SUMSTATS} ($(wc -l < "${SUMSTATS}") lines)"

# ---------------------------------------------------------------------------
# Standardise to GWAMA format
# ---------------------------------------------------------------------------
echo ""
echo "--- Cohort QC and GWAMA conversion ---"
if ! command -v "${RSCRIPT_BIN}" >/dev/null 2>&1 && [ ! -x "${RSCRIPT_BIN}" ]; then
  echo "ERROR: Rscript not found (RSCRIPT_BIN='${RSCRIPT_BIN}')" >&2
  exit 1
fi
echo "Using Rscript: ${RSCRIPT_BIN}"

"${RSCRIPT_BIN}" "${QC_SCRIPT}" \
  --input "${SUMSTATS}" \
  --cohort "${POPULATIONID}" \
  --build "${BUILD}" \
  --trait-type "${TRAITTYPE}" \
  --output-prefix "${POPULATIONID}" \
  --output-dir "${OUTDIR}" \
  "${COL_ARGS[@]}"

GWAMA_FILE="${OUTDIR}/${POPULATIONID}.GWAMA.txt"

if [ ! -s "${GWAMA_FILE}" ]; then
  echo "ERROR: QC step did not produce ${GWAMA_FILE}" >&2
  exit 1
fi

echo ""
echo "========================================"
echo "Site ${POPULATIONID} complete"
echo "========================================"
echo "GWAMA input : ${GWAMA_FILE} ($(( $(wc -l < "${GWAMA_FILE}") - 1 )) variants)"
echo "QC summary  : ${OUTDIR}/${POPULATIONID}.qc_summary.txt"
echo "QC log      : ${OUTDIR}/${POPULATIONID}.qc_log.txt"