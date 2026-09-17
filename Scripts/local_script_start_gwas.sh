

# Set the population id and program to use
# Values may come from the environment, or from the local JSON configuration file.
SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/test_configuration.json"

if [ -f "${CONFIG_FILE}" ]; then
  while IFS='=' read -r key value; do
    case "${key}" in
      PROGRAM)
        export PROGRAM="${value}"
        ;;
      POPULATIONID)
        export POPULATIONID="${value}"
        ;;
    esac
  done < <(python3 - "${CONFIG_FILE}" <<'PY'
import json, shlex, sys
cfg_path = sys.argv[1]
with open(cfg_path, encoding='utf-8') as f:
    cfg = json.load(f)
for key, env_name in {"program": "PROGRAM", "populationid": "POPULATIONID"}.items():
    value = cfg.get(key, "")
    if value is None:
        value = ""
    print(f"{env_name}={shlex.quote(str(value))}")
PY
)
fi

export POPULATIONID="${POPULATIONID:-}"
export PROGRAM="${PROGRAM:-}"

# site specific variables

#Use Data path from environment variable if set, otherwise use default
export DATA_PATH_DEF="Data/"


if [ -z "${DATAPATH:-}" ]; then
  echo "DATAPATH is not set in the environment; using default DATA_PATH=${DATA_PATH_DEF}"
  export DATAPATH="${DATA_PATH}"
else
  echo "DATAPATH is set to: ${DATAPATH}"
fi

export DATA_PATH


#Setting up the environment for GWAS analysis
export POPULATIONPATH="${DATA_PATH}/${POPULATIONID}/"

if [ "${PROGRAM}" = "regenie" ]; then
  echo "Running Regenie workflow"
  # Add Regenie-specific commands here
  echo "PROGRAM value: ${PROGRAM} is not yet supported"

elif [ "${PROGRAM}" = "plink" ] ; then
  echo "Running PLINK workflow"
  # Add PLINK-specific commands here
  echo "PROGRAM value: ${PROGRAM} is not yet supported"

elif [ "${PROGRAM}" = "gcta" ] ; then
  echo "Running GCTA workflow"
  # Add GCTA-specific commands here
  echo "PROGRAM value: ${PROGRAM} is not yet supported"

elif [ "${PROGRAM}" = "saige" ]; then
  echo "Running SAIGE workflow"
  echo "PROGRAM value: ${PROGRAM} is not yet supported"

  # Add SAIGE-specific commands here

else
  echo "Custom PROGRAM value: ${PROGRAM} is not yet supported"

  exit 1
fi



