#!/bin/bash
# One-shot setup + run: generates the 3 synthetic ancestry sites (EUR, EAS, AFR).
# No meta-analysis here -- just the genotype/phenotype/GWAS generation step.
#
# Usage:
#   1. scp requirements.txt, ancestry_sim.py, generate_federated_sites.sh and this
#      script into one directory on the server (they can sit flat, next to each other).
#   2. chmod +x run_on_server.sh && ./run_on_server.sh
set -euo pipefail
cd "$(dirname "$0")"

# ---- Parameters: the only place to change sample/SNP counts and causal SNPs ----
export SAMPLES="${SAMPLES:-100000 95000 110000}"    # per site, one per ancestry
export NSNPS="${NSNPS:-500000 480000 520000}"
export ANCESTRY="${ANCESTRY:-EUR EAS AFR}"
export SHARED_CAUSALS="${SHARED_CAUSALS:-20}"        # causal in every ancestry
export UNIQUE_CAUSALS="${UNIQUE_CAUSALS:-0}"         # causal in one ancestry only
export RUN_GWAS="${RUN_GWAS:-1}"
export PYTHON="python3"

# ---- 1. Lay out the expected directory structure ----
mkdir -p tools scripts
[ -f ancestry_sim.py ] && mv -f ancestry_sim.py tools/
[ -f generate_federated_sites.sh ] && mv -f generate_federated_sites.sh scripts/
chmod +x scripts/generate_federated_sites.sh

# ---- 2. Python virtual environment ----
if [ ! -d .venv ]; then
  python3 -m venv .venv
fi
source .venv/bin/activate
pip install -q -r requirements.txt

# ---- 3. LDAK binary, matched to this machine's OS ----
case "$(uname -s)" in
  Linux*)  ldak_asset="ldak6.3.linux"; ldak_name="ldak6.1.linux" ;;
  Darwin*) ldak_asset="ldak6.3.mac";   ldak_name="ldak6.1.mac" ;;
  *) echo "Unsupported OS: $(uname -s)"; exit 1 ;;
esac
if [ ! -f "tools/${ldak_name}" ]; then
  echo "Downloading LDAK (${ldak_asset})..."
  curl -sSL -o "tools/${ldak_name}" "https://github.com/dougspeed/LDAK/raw/main/${ldak_asset}"
  chmod +x "tools/${ldak_name}"
fi

# ---- 4. Sanity checks before a potentially long, large run ----
echo "SAMPLES=${SAMPLES}"; echo "NSNPS=${NSNPS}"; echo "ANCESTRY=${ANCESTRY}"
echo "cores: $(nproc 2>/dev/null || sysctl -n hw.ncpu)"
df -h . | tail -1

# ---- 5. Generate the sites (no meta-analysis) ----
n=$(echo "${SAMPLES}" | wc -w | tr -d ' ')
for s in $(seq 1 "${n}"); do
  ./scripts/generate_federated_sites.sh "${s}" || { echo "Site ${s} failed"; exit 1; }
done

echo "Done. Output in data/simulated_sites/site1 .. site${n}"
