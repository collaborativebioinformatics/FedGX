

# Set the population id and program to use
#TODO: these value comes from the serverrequest
export POPULATIONID
export PROGRAM

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

elif [ "${PROGRAM}" = "plink" ] ; then
  echo "Running PLINK workflow"
  # Add PLINK-specific commands here

elif [ "${PROGRAM}" = "gcta" ] ; then
  echo "Running GCTA workflow"
  # Add GCTA-specific commands here

elif [ "${PROGRAM}" = "saige" ]; then
  echo "Running SAIGE workflow"
  # Add SAIGE-specific commands here

elif [ "${PROGRAM}" = "custom" ] || [ "${PROGRAM}" = "CUSTOM" ]; then
  echo "Running custom workflow"
  # Add custom pipeline commands here

else
  echo "Unsupported PROGRAM value: ${PROGRAM}"
  exit 1
fi



