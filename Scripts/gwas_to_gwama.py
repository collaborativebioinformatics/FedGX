#!/usr/bin/env python3

# Convert REGENIE, SAIGE, or PLINK2 (--glm) GWAS output to GWAMA input format.
#
# Final effect allele (EA) is always the lower-frequency allele.
#
# MARKERNAME:
#     CHR:POS:EA:NEA
#
# ------------------------------------------------------------
# Usage
# ------------------------------------------------------------
#
# REGENIE case/control:
# python3 gwas_to_gwama.py input.regenie output.txt or regenie
#
# SAIGE case/control:
# python3 gwas_to_gwama.py input.saige output.txt or saige
#
# PLINK2 case/control:
# python3 gwas_to_gwama.py input.glm.logistic.hybrid output.txt or plink
#
# Quantitative trait:
# python3 gwas_to_gwama.py input.txt output.txt qt regenie
# python3 gwas_to_gwama.py input.txt output.txt qt saige
# python3 gwas_to_gwama.py input.glm.linear output.txt qt plink
#
#
# ------------------------------------------------------------
# Expected effect allele / frequency columns
# ------------------------------------------------------------
#
# REGENIE:
#   CHROM
#   GENPOS
#   ALLELE0
#   ALLELE1
#   A1FREQ
#   BETA
#   SE
#
#   Effect allele = ALLELE1
#   Effect allele frequency = A1FREQ
#
#
# SAIGE:
#   CHR
#   POS
#   Allele1
#   Allele2
#   AF_Allele2
#   BETA
#   SE
#
#   Effect allele = Allele2
#   Effect allele frequency = AF_Allele2
#
#
# PLINK2 --glm:
#   #CHROM / CHROM
#   POS
#   A1
#   A1_FREQ
#   OMITTED or REF/ALT1
#   BETA or OR
#   SE or LOG(OR)_SE
#
#   Effect allele = A1
#   Effect allele frequency = A1_FREQ
#
# ------------------------------------------------------------


import argparse
import pandas as pd
import numpy as np
from scipy.stats import norm


def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Convert REGENIE, SAIGE, or PLINK2 GWAS results "
            "to GWAMA input format."
        )
    )

    parser.add_argument(
        "input",
        help="Path to GWAS result file"
    )

    parser.add_argument(
        "output",
        help="Path to write GWAMA-formatted output"
    )

    parser.add_argument(
        "mode",
        choices=["or", "qt"],
        help=(
            "'or' = binary/case-control trait; "
            "'qt' = quantitative trait"
        )
    )

    parser.add_argument(
        "software",
        choices=["regenie", "saige", "plink"],
        help="GWAS software: regenie, saige, or plink"
    )

    return parser.parse_args()


def gwas2gwama(input_file, output_file, mode, software):

    # ------------------------------------------------------------
    # Read input
    # ------------------------------------------------------------

    gwas = pd.read_csv(
        input_file,
        sep=r"\s+"
    )


    # ------------------------------------------------------------
    # REGENIE
    # ------------------------------------------------------------

    if software == "regenie":

        required = [
            "CHROM",
            "GENPOS",
            "ALLELE0",
            "ALLELE1",
            "A1FREQ",
            "N",
            "BETA",
            "SE"
        ]

        missing = [
            col for col in required
            if col not in gwas.columns
        ]

        if missing:
            raise ValueError(
                "Missing REGENIE columns: "
                + ", ".join(missing)
            )

        gwas["CHR_STD"] = gwas["CHROM"]
        gwas["POS_STD"] = gwas["GENPOS"]

        # ALLELE1 is the original effect allele
        gwas["A0"] = gwas["ALLELE0"]
        gwas["A1"] = gwas["ALLELE1"]

        # Frequency of ALLELE1
        gwas["A1_FREQ"] = gwas["A1FREQ"]

        # Sample size
        gwas["N_STD"] = gwas["N"]

        # Effect estimate
        gwas["BETA_STD"] = gwas["BETA"]

        # Standard error
        gwas["SE_STD"] = gwas["SE"]


    # ------------------------------------------------------------
    # SAIGE
    # ------------------------------------------------------------

    elif software == "saige":

        required = [
            "CHR",
            "POS",
            "Allele1",
            "Allele2",
            "AF_Allele2",
            "BETA",
            "SE"
        ]

        missing = [
            col for col in required
            if col not in gwas.columns
        ]

        if missing:
            raise ValueError(
                "Missing SAIGE columns: "
                + ", ".join(missing)
            )

        gwas["CHR_STD"] = gwas["CHR"]
        gwas["POS_STD"] = gwas["POS"]

        # Allele2 is the original effect allele
        gwas["A0"] = gwas["Allele1"]
        gwas["A1"] = gwas["Allele2"]

        # Frequency of Allele2
        gwas["A1_FREQ"] = gwas["AF_Allele2"]

        # Sample size
        if "N" in gwas.columns:
            gwas["N_STD"] = gwas["N"]
        elif "N_case" in gwas.columns and "N_ctrl" in gwas.columns:
            gwas["N_STD"] = gwas["N_case"] + gwas["N_ctrl"]
        else:
            raise ValueError(
                "SAIGE input must contain N, or both N_case and N_ctrl, "
                "to write the GWAMA N column."
            )

        # Effect estimate
        gwas["BETA_STD"] = gwas["BETA"]

        # Standard error
        gwas["SE_STD"] = gwas["SE"]


    # ------------------------------------------------------------
    # PLINK2 --glm
    # ------------------------------------------------------------

    elif software == "plink":

        # PLINK2 normally uses #CHROM
        if "#CHROM" in gwas.columns:
            gwas["CHR_STD"] = gwas["#CHROM"]

        elif "CHROM" in gwas.columns:
            gwas["CHR_STD"] = gwas["CHROM"]

        else:
            raise ValueError(
                "PLINK input must contain #CHROM or CHROM."
            )


        # Position
        if "POS" not in gwas.columns:
            raise ValueError(
                "PLINK input must contain POS."
            )

        gwas["POS_STD"] = gwas["POS"]


        # --------------------------------------------------------
        # Keep additive GWAS result
        # --------------------------------------------------------

        if "TEST" in gwas.columns:

            n_before = len(gwas)

            gwas = gwas[
                gwas["TEST"] == "ADD"
            ].copy()

            print(
                f"PLINK ADD rows retained: "
                f"{len(gwas)} / {n_before}"
            )

            if len(gwas) == 0:
                raise ValueError(
                    "No TEST=ADD rows found in PLINK input."
                )


        # --------------------------------------------------------
        # Effect allele
        # --------------------------------------------------------

        if "A1" not in gwas.columns:
            raise ValueError(
                "PLINK input must contain A1."
            )

        gwas["A1"] = gwas["A1"]


        # --------------------------------------------------------
        # Effect allele frequency
        # --------------------------------------------------------

        if "A1_FREQ" not in gwas.columns:
            raise ValueError(
                "PLINK input must contain A1_FREQ. "
                "A1_FREQ is required to determine the "
                "lower-frequency allele."
            )

        gwas["A1_FREQ"] = gwas["A1_FREQ"]


        # --------------------------------------------------------
        # Determine the other allele
        # --------------------------------------------------------

        if "OMITTED" in gwas.columns:

            # Best PLINK2 representation of the non-A1 allele
            gwas["A0"] = gwas["OMITTED"]


        elif (
            "REF" in gwas.columns
            and "ALT1" in gwas.columns
        ):

            # For biallelic variants:
            #
            # If A1 == REF:
            #     other allele = ALT1
            #
            # Otherwise:
            #     other allele = REF

            gwas["A0"] = np.where(
                gwas["A1"] == gwas["REF"],
                gwas["ALT1"],
                gwas["REF"]
            )


        else:

            raise ValueError(
                "PLINK input must contain either OMITTED "
                "or both REF and ALT1."
            )


        # --------------------------------------------------------
        # Sample size
        # --------------------------------------------------------

        if "OBS_CT" in gwas.columns:
            gwas["N_STD"] = gwas["OBS_CT"]
        elif "N" in gwas.columns:
            gwas["N_STD"] = gwas["N"]
        else:
            raise ValueError(
                "PLINK input must contain OBS_CT or N "
                "to write the GWAMA N column."
            )


        # --------------------------------------------------------
        # PLINK effect estimate
        # --------------------------------------------------------

        if "BETA" in gwas.columns:

            # PLINK output already contains beta/log(OR)
            gwas["BETA_STD"] = gwas["BETA"]


        elif "OR" in gwas.columns:

            # For logistic regression:
            #
            # BETA = log(OR)

            if mode != "or":
                raise ValueError(
                    "PLINK OR column can only be used "
                    "with mode='or'."
                )

            if (gwas["OR"] <= 0).any():
                raise ValueError(
                    "PLINK OR values must be greater than 0."
                )

            gwas["BETA_STD"] = np.log(
                gwas["OR"]
            )


        else:

            raise ValueError(
                "PLINK input must contain BETA or OR."
            )


        # --------------------------------------------------------
        # PLINK standard error
        # --------------------------------------------------------

        if "SE" in gwas.columns:

            gwas["SE_STD"] = gwas["SE"]


        elif "LOG(OR)_SE" in gwas.columns:

            gwas["SE_STD"] = gwas["LOG(OR)_SE"]


        else:

            raise ValueError(
                "PLINK input must contain SE "
                "or LOG(OR)_SE."
            )


    # ------------------------------------------------------------
    # Convert numeric columns
    # ------------------------------------------------------------

    gwas["A1_FREQ"] = pd.to_numeric(
        gwas["A1_FREQ"],
        errors="coerce"
    )

    gwas["BETA_STD"] = pd.to_numeric(
        gwas["BETA_STD"],
        errors="coerce"
    )

    gwas["SE_STD"] = pd.to_numeric(
        gwas["SE_STD"],
        errors="coerce"
    )

    gwas["N_STD"] = pd.to_numeric(
        gwas["N_STD"],
        errors="coerce"
    )


    # ------------------------------------------------------------
    # Remove rows with missing required values
    # ------------------------------------------------------------

    n_before_missing = len(gwas)

    gwas = gwas.dropna(
        subset=[
            "CHR_STD",
            "POS_STD",
            "A0",
            "A1",
            "A1_FREQ",
            "BETA_STD",
            "SE_STD",
            "N_STD"
        ]
    ).copy()

    n_removed_missing = (
        n_before_missing - len(gwas)
    )


    # ------------------------------------------------------------
    # Check allele frequencies
    # ------------------------------------------------------------

    invalid_freq = (
        (gwas["A1_FREQ"] < 0)
        |
        (gwas["A1_FREQ"] > 1)
    )

    if invalid_freq.any():

        raise ValueError(
            f"{invalid_freq.sum()} variants have allele "
            "frequencies outside the range 0-1."
        )


    # ------------------------------------------------------------
    # Make lower-frequency allele the effect allele
    # ------------------------------------------------------------

    flip = gwas["A1_FREQ"] > 0.5


    # Original effect allele has frequency <= 0.5:
    #
    #     EA  = A1
    #     NEA = A0
    #
    # Original effect allele has frequency > 0.5:
    #
    #     EA  = A0
    #     NEA = A1

    gwas["EA"] = np.where(
        flip,
        gwas["A0"],
        gwas["A1"]
    )

    gwas["NEA"] = np.where(
        flip,
        gwas["A1"],
        gwas["A0"]
    )


    # ------------------------------------------------------------
    # Flip BETA when effect allele changes
    # ------------------------------------------------------------

    gwas["BETA"] = gwas["BETA_STD"]

    gwas.loc[
        flip,
        "BETA"
    ] = -gwas.loc[
        flip,
        "BETA"
    ]


    # ------------------------------------------------------------
    # Standard error does NOT change when allele is flipped
    # ------------------------------------------------------------

    gwas["SE"] = gwas["SE_STD"]


    # ------------------------------------------------------------
    # Final effect allele frequency
    # ------------------------------------------------------------

    gwas["EAF"] = np.where(
        flip,
        1 - gwas["A1_FREQ"],
        gwas["A1_FREQ"]
    )

    # Final sample size
    gwas["N"] = gwas["N_STD"]


    # ------------------------------------------------------------
    # Construct marker name
    #
    # CHR:POS:EA:NEA
    # ------------------------------------------------------------

    gwas["MARKERNAME"] = (
        gwas["CHR_STD"].astype(str)
        + ":"
        + gwas["POS_STD"].astype(str)
        + ":"
        + gwas["EA"].astype(str)
        + ":"
        + gwas["NEA"].astype(str)
    )


    # ------------------------------------------------------------
    # Binary / case-control trait
    # ------------------------------------------------------------

    if mode == "or":

        z_score = norm.ppf(0.975)

        # BETA is log(OR)
        gwas["OR"] = np.exp(
            gwas["BETA"]
        )

        gwas["OR_95L"] = np.exp(
            gwas["BETA"]
            - z_score * gwas["SE"]
        )

        gwas["OR_95U"] = np.exp(
            gwas["BETA"]
            + z_score * gwas["SE"]
        )


        output_columns = [
            "MARKERNAME",
            "EA",
            "NEA",
            "OR",
            "OR_95L",
            "OR_95U",
            "BETA",
            "SE",
            "EAF",
            "N"
        ]


    # ------------------------------------------------------------
    # Quantitative trait
    # ------------------------------------------------------------

    else:

        output_columns = [
            "MARKERNAME",
            "EA",
            "NEA",
            "BETA",
            "SE",
            "EAF",
            "N"
        ]


    # ------------------------------------------------------------
    # Write GWAMA output
    # ------------------------------------------------------------

    gwas[
        output_columns
    ].to_csv(
        output_file,
        sep="\t",
        index=False
    )


    # ------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------

    print()
    print("Conversion complete")
    print("-------------------")

    print(
        f"Software: {software}"
    )

    print(
        f"Trait mode: {mode}"
    )

    print(
        f"Variants written: {len(gwas)}"
    )

    print(
        f"Alleles flipped: {int(flip.sum())}"
    )

    print(
        f"Alleles not flipped: {int((~flip).sum())}"
    )

    print(
        f"Rows removed due to missing values: "
        f"{n_removed_missing}"
    )

    print(
        f"Output: {output_file}"
    )


def main():

    args = parse_args()

    gwas2gwama(
        args.input,
        args.output,
        args.mode,
        args.software
    )


if __name__ == "__main__":
    main()

