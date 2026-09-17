#!/bin/bash
# Generate a single federated learning site with genomic data for Parkinson's disease
# Usage: ./scripts/generate_federated_sites.sh <site_number>
# Example: ./scripts/generate_federated_sites.sh 1

# One site per ancestry: site 1 EUR, site 2 EAS, site 3 AFR
samples=(100000 95000 110000)
nsnps=(500000 480000 520000)
ancestry=(EUR EAS AFR)
nsites=${#samples[@]}
# Causal SNPs shared by all ancestries, and unique to each ancestry
# (one number for every ancestry, or per ancestry, e.g. UNIQUE_CAUSALS=EUR=5,EAS=5,AFR=10)
SHARED_CAUSALS="${SHARED_CAUSALS:-20}"
UNIQUE_CAUSALS="${UNIQUE_CAUSALS:-0}"

# Check if site number is provided
if [ $# -eq 0 ]; then
  echo "Error: Site number required"
  echo "Usage: ./scripts/generate_federated_sites.sh <site_number>"
  echo "Example: ./scripts/generate_federated_sites.sh 1"
  echo "Site number must be between 1 and ${nsites}"
  exit 1
fi

site=$1

# Validate site number
if ! [[ "$site" =~ ^[0-9]+$ ]] || [ "$site" -lt 1 ] || [ "$site" -gt ${nsites} ]; then
  echo "Error: Site number must be between 1 and ${nsites}"
  exit 1
fi

# LDAK binary path
case "$(uname -s)" in
    Darwin)
        export OS_TYPE="mac"
        ;;
    Linux*)
        export OS_TYPE="linux"
        ;;
    *)
        echo "Unsupported OS. Please use macOS or Linux."
        exit 1
        ;;
esac

LDAK="./tools/ldak6.1.${OS_TYPE}"
SIM="${PYTHON:-python3} ./tools/ancestry_sim.py"
REF_DIR="./data/simulated_sites/reference"

# Create output directory
OUTPUT_DIR="./data/simulated_sites/site${site}"
mkdir -p ${OUTPUT_DIR}

idx=$((site-1))
anc=${ancestry[$idx]}

echo "========================================"
echo "Federated Genomic Data Simulation"
echo "========================================"
echo ""
echo "Generating site ${site}: ${samples[$idx]} samples, ${nsnps[$idx]} SNPs, ancestry ${anc}"
echo "========================================"

# Shared SNPs, causal SNPs and per-ancestry effects (built once for all sites)
causal_settings="shared=${SHARED_CAUSALS} unique=${UNIQUE_CAUSALS}"
if [ ! -f ${REF_DIR}/reference.npz ]; then
  echo "Step 0: Building shared multi-ancestry reference (${causal_settings})..."
  $SIM reference --out ${REF_DIR} --num-shared-causals ${SHARED_CAUSALS} \
    --num-unique-causals ${UNIQUE_CAUSALS} --power -0.25 --num-phenos 1 || exit 1
  echo "${causal_settings}" > ${REF_DIR}/settings.txt
elif [ "$(cat ${REF_DIR}/settings.txt 2>/dev/null)" != "${causal_settings}" ]; then
  echo "Error: ${REF_DIR} was built with different causal settings; delete it to rebuild"
  exit 1
fi
num_causals=$(head -1 ${REF_DIR}/causals_${anc}.txt | wc -w | tr -d ' ')

# Generate ancestry-specific genotypes (also creates .covar file)
echo "Step 1: Generating ${anc} genotypes..."
$SIM site \
  --reference ${REF_DIR} \
  --ancestry ${anc} \
  --out ${OUTPUT_DIR}/site${site}_geno \
  --num-samples ${samples[$idx]} \
  --num-snps ${nsnps[$idx]} \
  --seed ${site}

# Generate Parkinson's disease phenotypes
echo "Step 2: Generating Parkinson's phenotypes..."
$LDAK \
  --make-phenos ${OUTPUT_DIR}/site${site}_pheno \
  --bfile ${OUTPUT_DIR}/site${site}_geno \
  --her 0.25 \
  --prevalence 0.01 \
  --num-causals ${num_causals} \
  --causals ${REF_DIR}/causals_${anc}.txt \
  --effects ${REF_DIR}/effects_${anc}.txt \
  --power -0.25 \
  --num-phenos 1 \
  --covar ${OUTPUT_DIR}/site${site}_geno.covar \
  --covar-her 0.1

# Optional per-site GWAS for meta-analysis
if [ "${RUN_GWAS:-0}" = "1" ]; then
  echo "Step 3: Per-site GWAS..."
  $SIM gwas --bfile ${OUTPUT_DIR}/site${site}_geno --pheno ${OUTPUT_DIR}/site${site}_pheno.pheno \
    --covar ${OUTPUT_DIR}/site${site}_geno.covar --causals ${REF_DIR}/causals_${anc}.txt --out ${OUTPUT_DIR}/site${site}.gwas.tsv
fi

echo ""
echo "========================================"
echo "Site ${site} generated successfully!"
echo "========================================"
echo ""
echo "Generated files:"
echo "  - Genotypes: ${OUTPUT_DIR}/site${site}_geno.bed/bim/fam"
echo "  - Phenotypes: ${OUTPUT_DIR}/site${site}_pheno.pheno"
echo "  - Covariates: ${OUTPUT_DIR}/site${site}_geno.covar"