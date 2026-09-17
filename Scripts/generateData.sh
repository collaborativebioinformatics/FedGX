#!/bin/bash
# Generate a single federated learning site with genomic data for Parkinson's disease.
#
# Three sites, one "ancestry" each, with a shared set of causal variants
# (so the meta-analysis has something to find) plus site-specific causal
# variants (so it also has heterogeneity).
#
# Usage:   ./scripts/generate_federated_sites.sh <site_number>
# Example: ./scripts/generate_federated_sites.sh 1

set -euo pipefail

# ---------------------------------------------------------------------------
# Site configuration
# ---------------------------------------------------------------------------
# Per-site variant counts. These differ, but they are NOT passed to --make-snps
# directly. LDAK spreads --num-snps evenly over 22 chromosomes, so calling it
# with 500000 at one site and 480000 at another puts the Nth variant on a
# different chromosome and at a different position at each site. Marker names
# would still line up, GWAMA would meta-analyse them without complaint, and the
# Manhattan plot would be meaningless.
#
# Instead every site generates the same MASTER grid (the largest count), which
# fixes a common chr/pos/name layout, and is then subsetted down to its own
# count. Different sites drop different variants, so coverage overlaps heavily
# but is not identical -- which is what real cohorts on different arrays or
# imputation panels look like.
samples=(100000 95000 110000)
nsnps=(500000 480000 520000)
ancestry=(EUR EAS AFR)

# Crude ancestry proxy: different MAF spectra per site. This is NOT real
# ancestry (see the caveats at the bottom) but it does give each site a
# different allele-frequency distribution, which is what the harmonisation code
# and GWAMA's allele-frequency check actually care about.
maf_low=(0.01 0.01 0.005)
maf_high=(0.50 0.45 0.50)

# Within-site population structure, so the auto-generated .covar file has
# something real to correct for. 1 = homogeneous.
populations=(1 2 3)

# Transferability of the shared effects to each ancestry. Site 1 is the
# reference; attenuation at the others is what produces non-zero I^2 in the
# GWAMA output. Set all to 1.00 for a homogeneous demo.
#
# IMPORTANT: LDAK rescales the noise term so the phenotype hits --her exactly,
# so scaling *every* effect by the same factor is completely absorbed and
# changes nothing. The scale is therefore applied only to the first half of the
# shared causals, which shifts their effects relative to the others and does
# survive the rescaling. Per-site --her does the rest.
effect_scale=(1.00 0.85 0.70)

# Per-site heritability. Differs slightly so the per-allele effect sizes are
# not identical across sites even for the unattenuated causals.
her=(0.25 0.22 0.20)

# Causal architecture
N_SHARED=12       # same variants at every site, effects scaled per ancestry
N_SPECIFIC=8      # different variants at each site

PREVALENCE=0.01
POWER=-0.25
COVAR_HER=0.1

# Set to 0 to let LDAK sample effect sizes itself. Shared causal *positions*
# stay shared either way, but their effect sizes would then differ per site.
USE_FIXED_EFFECTS=1

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
if [ $# -eq 0 ]; then
  echo "Error: Site number required"
  echo "Usage: ./scripts/generate_federated_sites.sh <site_number>"
  echo "Site number must be between 1 and ${#samples[@]}"
  exit 1
fi

site=$1

if ! [[ "$site" =~ ^[0-9]+$ ]] || [ "$site" -lt 1 ] || [ "$site" -gt "${#samples[@]}" ]; then
  echo "Error: Site number must be between 1 and ${#samples[@]}"
  exit 1
fi

idx=$((site - 1))

# ---------------------------------------------------------------------------
# LDAK binary
# ---------------------------------------------------------------------------
case "$(uname -s)" in
    Darwin) export OS_TYPE="mac" ;;
    Linux*) export OS_TYPE="linux" ;;
    *) echo "Unsupported OS. Please use macOS or Linux."; exit 1 ;;
esac

LDAK="./ldak6.3.${OS_TYPE}"

OUTPUT_DIR="./data/simulated_sites/site${site}"
mkdir -p "${OUTPUT_DIR}"

STEM="${OUTPUT_DIR}/site${site}_geno"
FULL="${STEM}_full"

# The master grid is the largest per-site count. Every site generates this many
# variants so that chromosome, position and name are identical everywhere.
MASTER_SNPS=0
for n in "${nsnps[@]}"; do
  [ "$n" -gt "$MASTER_SNPS" ] && MASTER_SNPS="$n"
done

echo "========================================"
echo "Federated Genomic Data Simulation"
echo "========================================"
echo "Site         : ${site}"
echo "Ancestry     : ${ancestry[$idx]}"
echo "Samples      : ${samples[$idx]}"
echo "SNPs         : ${nsnps[$idx]} of ${MASTER_SNPS} master grid"
echo "MAF range    : ${maf_low[$idx]} - ${maf_high[$idx]}"
echo "Populations  : ${populations[$idx]}"
echo "Causals      : ${N_SHARED} shared + ${N_SPECIFIC} site-specific"
echo "Effect scale : ${effect_scale[$idx]} (first half of shared causals)"
echo "Heritability : ${her[$idx]}"
echo "========================================"
echo ""

# ---------------------------------------------------------------------------
# Step 1: genotypes
# ---------------------------------------------------------------------------
echo "Step 1: Generating genotypes..."
$LDAK \
  --make-snps "${FULL}" \
  --num-samples "${samples[$idx]}" \
  --num-snps "${MASTER_SNPS}" \
  --maf-low "${maf_low[$idx]}" \
  --maf-high "${maf_high[$idx]}" \
  --populations "${populations[$idx]}"

# ---------------------------------------------------------------------------
# Step 2: pick the shared causal variants
# ---------------------------------------------------------------------------
# Chosen from the FULL bim at fixed fractional positions. Because every site
# generates the same number of variants over the same 22 chromosomes, this rule
# selects the same physical variants everywhere without any shared seed file or
# coordination step between site runs.
SHARED="${OUTPUT_DIR}/shared_causals.txt"

awk -v n="$N_SHARED" '
  { name[NR] = $2 }
  END { for (i = 1; i <= n; i++) print name[int(NR * i / (n + 1))] }
' "${FULL}.bim" > "${SHARED}"

echo ""
echo "Step 2: Shared causal variants selected (${N_SHARED})"

# ---------------------------------------------------------------------------
# Step 3: subset the master grid down to this site's variant count
# ---------------------------------------------------------------------------
# Selection is evenly spread across the genome (Bresenham-style) rather than a
# random sample or a head/tail truncation, so every chromosome stays covered.
# The per-site offset makes each site drop a different set of variants.
# Shared causals are always retained, and the final count is exactly nsnps.
TARGET="${nsnps[$idx]}"

if [ "${TARGET}" -lt "${MASTER_SNPS}" ]; then
  echo ""
  echo "Step 3: Subsetting ${MASTER_SNPS} -> ${TARGET} variants..."
  KEEPLIST="${OUTPUT_DIR}/site${site}_keep.snps"

  awk -v want="${TARGET}" -v off="${site}" '
    NR == FNR { shared[$1]; ns++; next }
    { if ($2 in shared) { print $2; next } cand[++nc] = $2 }
    END {
      need = want - ns
      if (need < 0) { print "ERROR: target smaller than shared causal count" > "/dev/stderr"; exit 1 }
      acc = (off * 7919) % nc
      for (i = 1; i <= nc; i++) {
        acc += need
        if (acc >= nc) { acc -= nc; print cand[i] }
      }
    }
  ' "${SHARED}" "${FULL}.bim" > "${KEEPLIST}"

  $LDAK \
    --make-bed "${STEM}" \
    --bfile "${FULL}" \
    --extract "${KEEPLIST}"

  cp "${FULL}.covar" "${STEM}.covar"
else
  echo ""
  echo "Step 3: Site uses the full master grid."
  for ext in bed bim fam covar; do
    cp "${FULL}.${ext}" "${STEM}.${ext}"
  done
fi

N_KEPT=$(wc -l < "${STEM}.bim")
echo "  Variants retained: ${N_KEPT} (target ${TARGET})"

# Sanity check: every shared causal must have survived subsetting.
MISSING=$(awk 'NR == FNR { keep[$2]; next } !($1 in keep)' "${STEM}.bim" "${SHARED}" | wc -l)
if [ "${MISSING}" -ne 0 ]; then
  echo "ERROR: ${MISSING} shared causal variant(s) missing from ${STEM}.bim"
  exit 1
fi

# ---------------------------------------------------------------------------
# Step 4: build the causals and effects files
# ---------------------------------------------------------------------------
# LDAK expects one row per phenotype and one column per causal predictor, so
# with --num-phenos 1 these are single whitespace-separated lines.
CAUSALS="${OUTPUT_DIR}/site${site}.causals"
EFFECTS="${OUTPUT_DIR}/site${site}.effects_in"

awk -v nspec="$N_SPECIFIC" -v nshared="$N_SHARED" -v site="$site" \
    -v scale="${effect_scale[$idx]}" \
    -v causals_out="$CAUSALS" -v effects_out="$EFFECTS" '
  function shared_effect(i) {
    # deterministic, alternating sign, identical at every site before scaling
    return ((i % 2) ? 1 : -1) * (0.60 + 0.04 * i)
  }
  # first file: the shared causal names, in order.
  # The ancestry scale is applied to the first half only -- a uniform scale
  # across all causals would be absorbed by the LDAK rescaling to --her.
  NR == FNR {
    nc++; cname[nc] = $1; shared[$1]
    ceff[nc] = shared_effect(nc) * ((nc <= nshared / 2) ? scale : 1.0)
    next
  }
  # second file: this site is bim, indexed by line
  { name[FNR] = $2; total = FNR }
  END {
    # site-specific causals: disjoint band, offset by site, skipping shared ones
    placed = 0
    j = int(total * 0.5) + site * 977
    while (placed < nspec) {
      j += 1301
      while (j > total) j -= total
      if (name[j] == "" || (name[j] in shared) || (name[j] in used)) continue
      used[name[j]]
      nc++; cname[nc] = name[j]; ceff[nc] = 0.70 + 0.03 * placed
      placed++
    }
    for (i = 1; i <= nc; i++) printf "%s%s", cname[i], (i < nc ? " " : "\n") > causals_out
    for (i = 1; i <= nc; i++) printf "%.6f%s", ceff[i], (i < nc ? " " : "\n") > effects_out
    printf "  Total causal variants: %d\n", nc
  }
' "${SHARED}" "${STEM}.bim"

echo ""
echo "Step 4: Causal files written"
echo "  ${CAUSALS}"
echo "  ${EFFECTS}"

# ---------------------------------------------------------------------------
# Step 5: phenotypes
# ---------------------------------------------------------------------------
echo ""
echo "Step 5: Generating Parkinson's phenotypes..."

PHENO_ARGS=(
  --make-phenos "${OUTPUT_DIR}/site${site}_pheno"
  --bfile "${STEM}"
  --her "${her[$idx]}"
  --prevalence "${PREVALENCE}"
  --power "${POWER}"
  --num-phenos 1
  --causals "${CAUSALS}"
  --covar "${STEM}.covar"
  --covar-her "${COVAR_HER}"
)

if [ "${USE_FIXED_EFFECTS}" -eq 1 ]; then
  PHENO_ARGS+=(--effects "${EFFECTS}")
fi

$LDAK "${PHENO_ARGS[@]}"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "========================================"
echo "Site ${site} (${ancestry[$idx]}) generated successfully!"
echo "========================================"
echo ""
echo "Generated files:"
echo "  - Genotypes:    ${STEM}.bed/bim/fam"
echo "  - Phenotypes:   ${OUTPUT_DIR}/site${site}_pheno.pheno"
echo "  - Liabilities:  ${OUTPUT_DIR}/site${site}_pheno.liab"
echo "  - True effects: ${OUTPUT_DIR}/site${site}_pheno.effects"
echo "  - Covariates:   ${STEM}.covar"
echo "  - Causal list:  ${CAUSALS}"
echo "  - Shared list:  ${SHARED}"
echo ""
echo "Shared causal variants (must be identical across all three sites):"
sed 's/^/  /' "${SHARED}"