#!/usr/bin/env Rscript

# ============================================================
# GWAS PIPELINE: NAMED COMMAND-LINE CONFIGURATION + COHORT QC
#
# Example:
# Rscript gwas_cohort_qc.R \
#   --input cohort.regenie.gz \
#   --cohort COHORT1 \
#   --build GRCh38 \
#   --trait-type quantitative \
#   --info-threshold 0.30 \
#   --min-n 30 \
#   --min-mac 6 \
#   --remove-palindromic TRUE \
#   --autosomes-only TRUE \
#   --snp-only TRUE \
#   --remove-duplicates TRUE \
#   --filter-test TRUE \
#   --test-value ADD \
#   --col-chr CHROM \
#   --col-pos GENPOS \
#   --col-id ID \
#   --col-ea ALLELE1 \
#   --col-nea ALLELE0 \
#   --col-eaf A1FREQ \
#   --col-beta BETA \
#   --col-se SE \
#   --col-n N \
#   --col-info INFO \
#   --col-log10p LOG10P \
#   --col-test TEST \
#   --col-chisq CHISQ \
#   --output-prefix COHORT1
# ============================================================
required_packages <- c(
    "optparse",
    "data.table",
    "R.utils"
)

missing_packages <- required_packages[
    !vapply(
        required_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]

if (length(missing_packages) > 0) {
    stop(
        "Missing required R packages: ",
        paste(missing_packages, collapse = ", "),
        ". Please install them before running this workflow."
    )
}

required_packages <- c(
    "optparse",
    "data.table",
    "R.utils"
)

missing_packages <- required_packages[
    !vapply(
        required_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]

if (length(missing_packages) > 0) {
    stop(
        "Missing required R packages: ",
        paste(missing_packages, collapse = ", "),
        ". Please install them before running this workflow."
    )
}

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(R.utils)
})
# ============================================================
# COMMAND-LINE ARGUMENTS
# ============================================================

option_list <- list(
    make_option("--input", type="character",
                help="Input GWAS summary-statistics file"),
    make_option("--cohort", type="character",
                help="Cohort name"),
    make_option("--build", type="character", default=NA_character_,
                help="Genome build, e.g. GRCh37 or GRCh38"),

    make_option("--trait-type", dest="trait_type", type="character", default="quantitative",
                help="Trait type for GWAMA output: quantitative or binary [default %default]"),

    make_option("--info-threshold", dest="info_threshold", type="double", default=NA_real_,
                help="Minimum imputation quality threshold"),
    make_option("--min-n", dest="min_n", type="double", default=30,
                help="Minimum per-variant sample size [default %default]"),
    make_option("--min-mac", dest="min_mac", type="double", default=6,
                help="Minimum minor allele count [default %default]"),

    make_option("--remove-palindromic", dest="remove_palindromic", type="character", default="TRUE",
                help="Remove A/T and C/G SNPs: TRUE/FALSE [default %default]"),
    make_option("--autosomes-only", dest="autosomes_only", type="character", default="TRUE",
                help="Keep chromosomes 1-22 only: TRUE/FALSE [default %default]"),
    make_option("--snp-only", dest="snp_only", type="character", default="TRUE",
                help="Keep SNPs only: TRUE/FALSE [default %default]"),
    make_option("--remove-duplicates", dest="remove_duplicates", type="character", default="TRUE",
                help="Remove duplicate variants: TRUE/FALSE [default %default]"),
    make_option("--filter-test", dest="filter_test", type="character", default="FALSE",
                help="Filter association test type: TRUE/FALSE [default %default]"),
    make_option("--test-value", dest="test_value", type="character", default="ADD",
                help="Association test to retain [default %default]"),

    make_option("--col-chr", dest="col_chr", type="character",
                help="Chromosome column"),
    make_option("--col-pos", dest="col_pos", type="character",
                help="Position column"),
    make_option("--col-id", dest="col_id", type="character",
                help="Variant ID column"),
    make_option("--col-ea", dest="col_ea", type="character",
                help="Effect allele column"),
    make_option("--col-nea", dest="col_nea", type="character",
                help="Non-effect allele column"),
    make_option("--col-eaf", dest="col_eaf", type="character",
                help="Effect allele frequency column"),
    make_option("--col-beta", dest="col_beta", type="character",
                help="Beta/effect-size column"),
    make_option("--col-se", dest="col_se", type="character",
                help="Standard error column"),
    make_option("--col-n", dest="col_n", type="character",
                help="Sample-size column"),

    make_option("--col-info", dest="col_info", type="character", default=NA_character_,
                help="Imputation quality column; use NULL if unavailable"),
    make_option("--col-p", dest="col_p", type="character", default=NA_character_,
                help="P-value column; use NULL if unavailable"),
    make_option("--col-log10p", dest="col_log10p", type="character", default=NA_character_,
                help="-log10(P) column; use NULL if unavailable"),
    make_option("--col-test", dest="col_test", type="character", default=NA_character_,
                help="Association test column; use NULL if unavailable"),
    make_option("--col-chisq", dest="col_chisq", type="character", default=NA_character_,
                help="Chi-square statistic column; use NULL if unavailable"),

    make_option("--output-prefix", dest="output_prefix", type="character", default=NA_character_,
                help="Output prefix; defaults to cohort name"),
    make_option("--output-dir", dest="output_dir", type="character", default=".",
                help="Output directory [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ============================================================
# ARGUMENT HELPERS / VALIDATION
# ============================================================

parse_bool <- function(x, name) {
    x <- toupper(trimws(x))
    if (x %in% c("TRUE", "T", "1", "YES", "Y")) return(TRUE)
    if (x %in% c("FALSE", "F", "0", "NO", "N")) return(FALSE)
    stop(name, " must be TRUE or FALSE. Received: ", x)
}

na_to_null <- function(x) {
    if (length(x) == 0 || is.na(x) || trimws(x) == "" ||
        toupper(trimws(x)) == "NULL") {
        return(NULL)
    }
    x
}

required_args <- c(
    input="input", cohort="cohort",
    col_chr="col_chr", col_pos="col_pos", col_id="col_id",
    col_ea="col_ea", col_nea="col_nea", col_eaf="col_eaf",
    col_beta="col_beta", col_se="col_se", col_n="col_n"
)

for (nm in names(required_args)) {
    value <- opt[[required_args[[nm]]]]
    if (is.null(value) || is.na(value) || value == "") {
        stop("Missing required argument: --", gsub("_", "-", required_args[[nm]]))
    }
}

remove_palindromic <- parse_bool(opt$remove_palindromic, "--remove-palindromic")
autosomes_only <- parse_bool(opt$autosomes_only, "--autosomes-only")
snp_only <- parse_bool(opt$snp_only, "--snp-only")
remove_duplicates <- parse_bool(opt$remove_duplicates, "--remove-duplicates")
filter_test <- parse_bool(opt$filter_test, "--filter-test")

trait_type <- tolower(trimws(opt$trait_type))
if (!trait_type %in% c("quantitative", "binary")) {
    stop("--trait-type must be either 'quantitative' or 'binary'. Received: ", opt$trait_type)
}

col_info   <- na_to_null(opt$col_info)
col_p      <- na_to_null(opt$col_p)
col_log10p <- na_to_null(opt$col_log10p)
col_test   <- na_to_null(opt$col_test)
col_chisq  <- na_to_null(opt$col_chisq)

if (is.null(col_p) && is.null(col_log10p)) {
    stop("Specify either --col-p or --col-log10p.")
}

if (!is.null(col_info) && is.na(opt$info_threshold)) {
    stop("--info-threshold is required when --col-info is supplied.")
}

output_prefix <- if (
    is.na(opt$output_prefix) || opt$output_prefix == ""
) opt$cohort else opt$output_prefix

dir.create(opt$output_dir, recursive=TRUE, showWarnings=FALSE)

# ============================================================
# QC FUNCTION
# ============================================================

run_gwas_qc <- function(
    data,
    cohort_name,
    genome_build = NA_character_,
    info_threshold = NULL,
    min_n = 30,
    min_mac = 6,
    autosomes_only = TRUE,
    remove_palindromic = TRUE,
    snp_only = TRUE,
    remove_duplicates = TRUE,
    filter_test = FALSE,
    test_value = "ADD",
    col_chr,
    col_pos,
    col_id,
    col_ea,
    col_nea,
    col_eaf,
    col_beta,
    col_se,
    col_n,
    col_info = NULL,
    col_p = NULL,
    col_log10p = NULL,
    col_test = NULL,
    col_chisq = NULL
) {

    # --------------------------------------------------------
    # Internal helpers
    # --------------------------------------------------------

    add_qc <- function(step, before, after) {
        data.table(
            QC_STEP = step,
            N_BEFORE = before,
            N_REMOVED = before - after,
            N_REMAINING = after,
            PERCENT_REMOVED = ifelse(
                before > 0,
                100 * (before - after) / before,
                NA_real_
            )
        )
    }

    check_column <- function(dat, column_name, label, required = TRUE) {

        if (is.null(column_name)) {

            if (required) {
                stop(label, " column mapping is NULL.")
            }

            return(invisible(FALSE))
        }

        if (!column_name %in% names(dat)) {

            if (required) {
                stop(
                    label,
                    " column not found: ",
                    column_name
                )
            }

            return(invisible(FALSE))
        }

        invisible(TRUE)
    }


    # --------------------------------------------------------
    # Prepare input
    # --------------------------------------------------------

    d0 <- as.data.table(data)

    if (nrow(d0) == 0) {
        stop("Input data contains zero rows.")
    }

    n_original <- nrow(d0)


    # --------------------------------------------------------
    # Validate required column mappings
    # --------------------------------------------------------

    check_column(d0, col_chr,  "Chromosome")
    check_column(d0, col_pos,  "Position")
    check_column(d0, col_id,   "Variant ID")
    check_column(d0, col_ea,   "Effect allele")
    check_column(d0, col_nea,  "Non-effect allele")
    check_column(d0, col_eaf,  "Effect allele frequency")
    check_column(d0, col_beta, "Beta")
    check_column(d0, col_se,   "SE")
    check_column(d0, col_n,    "N")

    if (is.null(col_p) && is.null(col_log10p)) {
        stop("Provide either col_p or col_log10p.")
    }

    if (!is.null(col_p)) {
        check_column(d0, col_p, "P-value")
    }

    if (!is.null(col_log10p)) {
        check_column(d0, col_log10p, "-log10(P)")
    }

    if (!is.null(col_info)) {
        check_column(d0, col_info, "Imputation quality")
    }

    if (filter_test) {
        if (is.null(col_test)) {
            stop("--filter-test TRUE requires --col-test.")
        }

        check_column(
            d0,
            col_test,
            "Test type",
            required = TRUE
        )
    } else if (!is.null(col_test)) {
        check_column(d0, col_test, "Test type", required = FALSE)
    }

    if (!is.null(col_chisq)) {
        check_column(d0, col_chisq, "Chi-square", required = FALSE)
    }


    # --------------------------------------------------------
    # Standardize core variables
    # --------------------------------------------------------

    d <- data.table(
        CHR = as.character(d0[[col_chr]]),

        POS = suppressWarnings(
            as.integer(d0[[col_pos]])
        ),

        SNP = as.character(d0[[col_id]]),

        EA = toupper(
            as.character(d0[[col_ea]])
        ),

        NEA = toupper(
            as.character(d0[[col_nea]])
        ),

        EAF = suppressWarnings(
            as.numeric(d0[[col_eaf]])
        ),

        BETA = suppressWarnings(
            as.numeric(d0[[col_beta]])
        ),

        SE = suppressWarnings(
            as.numeric(d0[[col_se]])
        ),

        N = suppressWarnings(
            as.numeric(d0[[col_n]])
        )
    )


    # --------------------------------------------------------
    # Optional INFO/imputation-quality field
    # --------------------------------------------------------

    if (!is.null(col_info)) {

        d[, INFO := suppressWarnings(
            as.numeric(d0[[col_info]])
        )]

    } else {

        d[, INFO := NA_real_]
    }


    # --------------------------------------------------------
    # Optional TEST field
    # --------------------------------------------------------

    if (
        !is.null(col_test) &&
        col_test %in% names(d0)
    ) {

        d[, TEST := as.character(d0[[col_test]])]

    } else {

        d[, TEST := NA_character_]
    }


    # --------------------------------------------------------
    # Optional CHISQ field
    # --------------------------------------------------------

    if (
        !is.null(col_chisq) &&
        col_chisq %in% names(d0)
    ) {

        d[, CHISQ := suppressWarnings(
            as.numeric(d0[[col_chisq]])
        )]

    } else {

        d[, CHISQ := NA_real_]
    }


    # --------------------------------------------------------
    # P-value handling
    # --------------------------------------------------------

    if (!is.null(col_p)) {

        d[, P := suppressWarnings(
            as.numeric(d0[[col_p]])
        )]

        # Mark truly invalid P-values as missing.
        d[P < 0 | P > 1, P := NA_real_]

        # Avoid Inf when software reports extremely small P-values as 0.
        # The reported P remains 0, while LOG10P is capped at machine precision.
        d[, LOG10P := -log10(
            pmax(P, .Machine$double.xmin)
        )]

    } else {

        d[, LOG10P := suppressWarnings(
            as.numeric(d0[[col_log10p]])
        )]

        # May numerically underflow for extremely small P.
        # LOG10P is retained and used for P-Z diagnostics.
        d[, P := 10^(-LOG10P)]
    }


    # --------------------------------------------------------
    # Initialize QC log
    # --------------------------------------------------------

    qc_log <- list()


    # --------------------------------------------------------
    # Filter association test type, if requested
    # --------------------------------------------------------

    if (
        filter_test &&
        !all(is.na(d$TEST))
    ) {

        before <- nrow(d)

        d <- d[
            !is.na(TEST) &
            TEST == test_value
        ]

        qc_log[[length(qc_log) + 1]] <- add_qc(
            paste0("Test != ", test_value),
            before,
            nrow(d)
        )
    }


    # --------------------------------------------------------
    # Remove missing essential values
    # --------------------------------------------------------

    before <- nrow(d)

    d <- d[
        !is.na(CHR) &
        !is.na(POS) &
        !is.na(SNP) &
        !is.na(EA) &
        !is.na(NEA) &
        !is.na(EAF) &
        !is.na(BETA) &
        !is.na(SE) &
        !is.na(N) &
        !is.na(LOG10P)
    ]

    qc_log[[length(qc_log) + 1]] <- add_qc(
        "Missing essential values",
        before,
        nrow(d)
    )


    # --------------------------------------------------------
    # Normalize chromosome labels
    # --------------------------------------------------------

    d[, CHR := gsub(
        "^chr",
        "",
        CHR,
        ignore.case = TRUE
    )]


    # --------------------------------------------------------
    # Remove invalid numeric values
    # --------------------------------------------------------

    before <- nrow(d)

    d <- d[
        is.finite(POS) &
        POS > 0 &
        is.finite(EAF) &
        EAF >= 0 &
        EAF <= 1 &
        is.finite(BETA) &
        is.finite(SE) &
        SE > 0 &
        is.finite(N) &
        N > 0 &
        is.finite(LOG10P) &
        LOG10P >= 0
    ]

    qc_log[[length(qc_log) + 1]] <- add_qc(
        "Invalid numeric values",
        before,
        nrow(d)
    )


    # --------------------------------------------------------
    # SNP-only allele filter
    # --------------------------------------------------------

    if (snp_only) {

        valid_alleles <- c("A", "C", "G", "T")

        before <- nrow(d)

        d <- d[
            EA %in% valid_alleles &
            NEA %in% valid_alleles &
            EA != NEA
        ]

        qc_log[[length(qc_log) + 1]] <- add_qc(
            "Invalid/non-SNP alleles",
            before,
            nrow(d)
        )
    }


    # --------------------------------------------------------
    # Remove monomorphic variants
    # --------------------------------------------------------

    before <- nrow(d)

    d <- d[
        EAF > 0 &
        EAF < 1
    ]

    qc_log[[length(qc_log) + 1]] <- add_qc(
        "Monomorphic variants",
        before,
        nrow(d)
    )


    # --------------------------------------------------------
    # Calculate MAF
    # --------------------------------------------------------

    d[, MAF := pmin(
        EAF,
        1 - EAF
    )]


    # --------------------------------------------------------
    # Sample-size filter
    # --------------------------------------------------------

    before <- nrow(d)

    d <- d[
        N >= min_n
    ]

    qc_log[[length(qc_log) + 1]] <- add_qc(
        paste0("N < ", min_n),
        before,
        nrow(d)
    )


    # --------------------------------------------------------
    # Approximate minor allele count
    #
    # Diploid autosomal approximation:
    # MAC = 2 * N * MAF
    # --------------------------------------------------------

    d[, MAC := 2 * N * MAF]

    before <- nrow(d)

    d <- d[
        MAC > min_mac
    ]

    qc_log[[length(qc_log) + 1]] <- add_qc(
        paste0("MAC <= ", min_mac),
        before,
        nrow(d)
    )


    # --------------------------------------------------------
    # Imputation-quality filter
    # --------------------------------------------------------

    if (!is.null(col_info)) {

        if (is.null(info_threshold)) {
            stop(
                "col_info was supplied but info_threshold is NULL."
            )
        }

        before <- nrow(d)

        d <- d[
            !is.na(INFO) &
            is.finite(INFO)
        ]

        qc_log[[length(qc_log) + 1]] <- add_qc(
            "Missing/non-finite INFO",
            before,
            nrow(d)
        )

        before <- nrow(d)

        d <- d[
            INFO >= info_threshold
        ]

        qc_log[[length(qc_log) + 1]] <- add_qc(
            paste0("INFO < ", info_threshold),
            before,
            nrow(d)
        )
    }


    # --------------------------------------------------------
    # Chromosome numeric representation
    # --------------------------------------------------------

    d[, CHR_NUM := suppressWarnings(
        as.integer(CHR)
    )]


    # --------------------------------------------------------
    # Autosomal filter
    # --------------------------------------------------------

    if (autosomes_only) {

        before <- nrow(d)

        d <- d[
            !is.na(CHR_NUM) &
            CHR_NUM >= 1 &
            CHR_NUM <= 22
        ]

        qc_log[[length(qc_log) + 1]] <- add_qc(
            "Non-autosomal variants",
            before,
            nrow(d)
        )
    }


    # --------------------------------------------------------
    # Flag palindromic SNPs
    #
    # A/T, T/A, C/G, G/C
    # --------------------------------------------------------

    d[, PALINDROMIC :=
        (EA == "A" & NEA == "T") |
        (EA == "T" & NEA == "A") |
        (EA == "C" & NEA == "G") |
        (EA == "G" & NEA == "C")
    ]


    # --------------------------------------------------------
    # Optional palindromic-SNP removal
    # --------------------------------------------------------

    if (remove_palindromic) {

        before <- nrow(d)

        d <- d[
            PALINDROMIC == FALSE
        ]

        qc_log[[length(qc_log) + 1]] <- add_qc(
            "Palindromic SNPs",
            before,
            nrow(d)
        )
    }


    # --------------------------------------------------------
    # Standardized variant identifiers
    # --------------------------------------------------------

    d[, CHR_POS := paste0(
        CHR,
        ":",
        POS
    )]

    # The key must not depend on which allele the tool happened to treat as
    # the effect allele. If one cohort reports EA=A/NEA=G and another
    # EA=G/NEA=A for the same variant, an EA:NEA key produces two different
    # strings, GWAMA sees two different markers, and every variant comes back
    # with n_studies = 1. Sorting the alleles in the key fixes that; the EA
    # and NEA columns still carry the direction, which is what GWAMA aligns on.
    d[, VARIANT_ID := paste(
        CHR,
        POS,
        pmin(EA, NEA),
        pmax(EA, NEA),
        sep = ":"
    )]


    # --------------------------------------------------------
    # Duplicate removal
    #
    # Keeps the row with:
    #   1. largest N
    #   2. largest INFO, if available
    # --------------------------------------------------------

    if (remove_duplicates) {

        before <- nrow(d)

        if (
            "INFO" %in% names(d) &&
            !all(is.na(d$INFO))
        ) {

            setorder(
                d,
                VARIANT_ID,
                -N,
                -INFO
            )

        } else {

            setorder(
                d,
                VARIANT_ID,
                -N
            )
        }

        d <- d[
            !duplicated(VARIANT_ID)
        ]

        qc_log[[length(qc_log) + 1]] <- add_qc(
            "Duplicate variants",
            before,
            nrow(d)
        )
    }


    # --------------------------------------------------------
    # Stop cleanly if no variants remain
    # --------------------------------------------------------

    if (nrow(d) == 0) {

        qc_table <- rbindlist(
            qc_log,
            fill = TRUE
        )

        qc_table[, COHORT := cohort_name]

        setcolorder(
            qc_table,
            c(
                "COHORT",
                "QC_STEP",
                "N_BEFORE",
                "N_REMOVED",
                "N_REMAINING",
                "PERCENT_REMOVED"
            )
        )

        return(
            list(
                clean_data = d,
                qc_summary = data.table(
                    COHORT = cohort_name,
                    GENOME_BUILD = genome_build,
                    N_ORIGINAL = n_original,
                    N_CLEAN = 0,
                    INFO_THRESHOLD = ifelse(
                        is.null(col_info),
                        NA_real_,
                        info_threshold
                    ),
                    MIN_N = min_n,
                    MIN_MAC = min_mac,
                    MEDIAN_N = NA_real_,
                    MAX_N = NA_real_,
                    MEDIAN_SE = NA_real_,
                    LAMBDA_GC = NA_real_,
                    MEDIAN_ABS_PZ_DIFF = NA_real_
                ),
                qc_log = qc_table,
                pz_summary = data.table(),
                se_n_summary = data.table()
            )
        )
    }


    # --------------------------------------------------------
    # P-Z consistency
    #
    # For SPA/Firth-corrected binary-trait tests (e.g. SAIGE),
    # disagreement between reported P and BETA/SE-derived P can
    # be expected and should not automatically be treated as QC
    # failure.
    # --------------------------------------------------------

    d[, Z := BETA / SE]

    d[, LOG_P_FROM_Z :=
        log(2) +
        pnorm(
            -abs(Z),
            log.p = TRUE
        )
    ]

    d[, LOG10P_FROM_Z :=
        -LOG_P_FROM_Z / log(10)
    ]

    d[, PZ_DIFF :=
        LOG10P - LOG10P_FROM_Z
    ]


    # --------------------------------------------------------
    # P-Z summary
    # --------------------------------------------------------

    pz_summary <- data.table(
        COHORT = cohort_name,
        N_VARIANTS = nrow(d),

        MEDIAN_ABS_PZ_DIFF = median(
            abs(d$PZ_DIFF),
            na.rm = TRUE
        ),

        P95_ABS_PZ_DIFF = as.numeric(
            quantile(
                abs(d$PZ_DIFF),
                probs = 0.95,
                na.rm = TRUE
            )
        ),

        MAX_ABS_PZ_DIFF = max(
            abs(d$PZ_DIFF),
            na.rm = TRUE
        ),

        CORRELATION = cor(
            d$LOG10P,
            d$LOG10P_FROM_Z,
            use = "complete.obs"
        )
    )


    # --------------------------------------------------------
    # Lambda GC
    #
    # Calculate from the reported association P-values so that
    # the diagnostic reflects the actual test used by the GWAS
    # software (important for SPA/Firth-corrected binary tests).
    #
    # For P=0 due to numerical underflow, use machine precision.
    # --------------------------------------------------------

    p_lambda <- d[
        !is.na(P) &
        is.finite(P) &
        P >= 0 &
        P <= 1,
        P
    ]

    p_lambda <- pmax(
        p_lambda,
        .Machine$double.xmin
    )

    chisq_values <- qchisq(
        p_lambda,
        df = 1,
        lower.tail = FALSE
    )

    lambda_gc <- median(
        chisq_values,
        na.rm = TRUE
    ) / qchisq(
        0.5,
        df = 1
    )


    # --------------------------------------------------------
    # SE-N summary
    # --------------------------------------------------------

    n_max <- max(
        d$N,
        na.rm = TRUE
    )

    n_median <- median(
        d$N,
        na.rm = TRUE
    )

    se_median <- median(
        d$SE,
        na.rm = TRUE
    )

    se_n_summary <- data.table(
        COHORT = cohort_name,
        N_MAX = n_max,
        N_MEDIAN = n_median,
        SQRT_N_MAX = sqrt(n_max),
        MEDIAN_SE = se_median,
        INVERSE_MEDIAN_SE = 1 / se_median
    )


    # --------------------------------------------------------
    # QC log table
    # --------------------------------------------------------

    qc_table <- rbindlist(
        qc_log,
        fill = TRUE
    )

    qc_table[, COHORT := cohort_name]

    setcolorder(
        qc_table,
        c(
            "COHORT",
            "QC_STEP",
            "N_BEFORE",
            "N_REMOVED",
            "N_REMAINING",
            "PERCENT_REMOVED"
        )
    )


    # --------------------------------------------------------
    # Overall cohort summary
    # --------------------------------------------------------

    qc_summary <- data.table(
        COHORT = cohort_name,
        GENOME_BUILD = genome_build,

        N_ORIGINAL = n_original,
        N_CLEAN = nrow(d),

        INFO_THRESHOLD = ifelse(
            is.null(col_info),
            NA_real_,
            info_threshold
        ),

        MIN_N = min_n,
        MIN_MAC = min_mac,

        MIN_EAF = min(
            d$EAF,
            na.rm = TRUE
        ),

        MEDIAN_EAF = median(
            d$EAF,
            na.rm = TRUE
        ),

        MIN_MAF = min(
            d$MAF,
            na.rm = TRUE
        ),

        MEDIAN_MAF = median(
            d$MAF,
            na.rm = TRUE
        ),

        MEDIAN_N = n_median,
        MAX_N = n_max,

        MEDIAN_SE = se_median,

        LAMBDA_GC = lambda_gc,

        MEDIAN_ABS_PZ_DIFF = median(
            abs(d$PZ_DIFF),
            na.rm = TRUE
        )
    )

    if (
        "INFO" %in% names(d) &&
        !all(is.na(d$INFO))
    ) {

        qc_summary[, `:=`(
            MIN_INFO = min(
                d$INFO,
                na.rm = TRUE
            ),

            MEDIAN_INFO = median(
                d$INFO,
                na.rm = TRUE
            )
        )]
    }


    # --------------------------------------------------------
    # Sort cleaned data
    # --------------------------------------------------------

    if (autosomes_only) {

        setorder(
            d,
            CHR_NUM,
            POS
        )

    } else {

        setorder(
            d,
            CHR,
            POS
        )
    }


    # --------------------------------------------------------
    # Return all objects
    # --------------------------------------------------------

    return(
        list(
            clean_data = d,
            qc_summary = qc_summary,
            qc_log = qc_table,
            pz_summary = pz_summary,
            se_n_summary = se_n_summary
        )
    )
}


# ============================================================
# OPTIONAL HELPER: P-Z PLOT
#
# Usage:
#   plot_pz(qc_result$clean_data, "COHORT1")
# ============================================================

plot_pz <- function(
    clean_data,
    cohort_name = ""
) {

    d <- as.data.table(clean_data)

    if (
        !all(
            c(
                "LOG10P",
                "LOG10P_FROM_Z"
            ) %in% names(d)
        )
    ) {
        stop(
            "clean_data must contain LOG10P and LOG10P_FROM_Z."
        )
    }

    plot(
        d$LOG10P_FROM_Z,
        d$LOG10P,
        pch = 16,
        cex = 0.25,

        xlab = expression(
            -log[10](P[beta/SE])
        ),

        ylab = expression(
            -log[10](P[reported])
        ),

        main = paste0(
            cohort_name,
            ": P-Z consistency"
        )
    )

    abline(
        a = 0,
        b = 1
    )
}


# ============================================================
# OPTIONAL HELPER: QQ PLOT
#
# Usage:
#   plot_qq(qc_result$clean_data, "COHORT1")
# ============================================================

plot_qq <- function(
    clean_data,
    cohort_name = ""
) {

    d <- as.data.table(clean_data)

    if (!"P" %in% names(d)) {
        stop("clean_data must contain P.")
    }

    qq_p <- d[
        is.finite(P) &
        P > 0 &
        P <= 1,
        P
    ]

    if (length(qq_p) == 0) {
        warning(
            "No numerically representable P-values available for QQ plot."
        )

        return(
            invisible(NULL)
        )
    }

    qq_p <- sort(qq_p)

    observed <- -log10(
        qq_p
    )

    expected <- -log10(
        ppoints(
            length(qq_p)
        )
    )

    plot(
        expected,
        observed,
        pch = 16,
        cex = 0.30,

        xlab = expression(
            Expected~~-log[10](P)
        ),

        ylab = expression(
            Observed~~-log[10](P)
        ),

        main = paste0(
            cohort_name,
            ": QQ plot"
        )
    )

    abline(
        a = 0,
        b = 1
    )
}




# ============================================================
# MAIN PIPELINE
# ============================================================

cat("\n========================================\n")
cat("GWAS QC CONFIGURATION\n")
cat("========================================\n")
cat("Cohort              :", opt$cohort, "\n")
cat("Input               :", opt$input, "\n")
cat("Genome build        :", opt$build, "\n")
cat("Trait type          :", trait_type, "\n")
cat("INFO threshold      :", ifelse(is.null(col_info), "not used", opt$info_threshold), "\n")
cat("Minimum N           :", opt$min_n, "\n")
cat("Minimum MAC         :", opt$min_mac, "\n")
cat("Remove palindromic  :", remove_palindromic, "\n")
cat("Autosomes only      :", autosomes_only, "\n")
cat("SNP only            :", snp_only, "\n")
cat("Remove duplicates   :", remove_duplicates, "\n")
cat("Filter test         :", filter_test, "\n")
cat("Test value          :", opt$test_value, "\n")
cat("Output directory    :", opt$output_dir, "\n")
cat("Output prefix       :", output_prefix, "\n")
cat("========================================\n\n")

cat("Reading input data...\n")
dat <- fread(opt$input)
cat("Rows read:", format(nrow(dat), big.mark=","), "\n")

# ------------------------------------------------------------
# INSERT ANY PRE-QC PIPELINE CODE HERE, IF NEEDED
# ------------------------------------------------------------

qc <- run_gwas_qc(
    data = dat,
    cohort_name = opt$cohort,
    genome_build = opt$build,

    info_threshold = if (is.null(col_info)) NULL else opt$info_threshold,
    min_n = opt$min_n,
    min_mac = opt$min_mac,

    autosomes_only = autosomes_only,
    remove_palindromic = remove_palindromic,
    snp_only = snp_only,
    remove_duplicates = remove_duplicates,

    filter_test = filter_test,
    test_value = opt$test_value,

    col_chr = opt$col_chr,
    col_pos = opt$col_pos,
    col_id = opt$col_id,

    col_ea = opt$col_ea,
    col_nea = opt$col_nea,
    col_eaf = opt$col_eaf,

    col_beta = opt$col_beta,
    col_se = opt$col_se,
    col_n = opt$col_n,

    col_info = col_info,
    col_p = col_p,
    col_log10p = col_log10p,
    col_test = col_test,
    col_chisq = col_chisq
)

clean <- qc$clean_data

# ============================================================
# WRITE QC OUTPUTS
# ============================================================

clean_file <- file.path(
    opt$output_dir,
    paste0(output_prefix, ".cleaned.txt")
)
summary_file <- file.path(
    opt$output_dir,
    paste0(output_prefix, ".qc_summary.txt")
)
log_file <- file.path(
    opt$output_dir,
    paste0(output_prefix, ".qc_log.txt")
)
pz_summary_file <- file.path(
    opt$output_dir,
    paste0(output_prefix, ".pz_summary.txt")
)
sen_file <- file.path(
    opt$output_dir,
    paste0(output_prefix, ".se_n_summary.txt")
)

fwrite(qc$clean_data, clean_file, sep="\t", quote=FALSE, na="NA")
fwrite(qc$qc_summary, summary_file, sep="\t", quote=FALSE, na="NA")
fwrite(qc$qc_log, log_file, sep="\t", quote=FALSE, na="NA")

if (nrow(qc$pz_summary) > 0) {
    fwrite(qc$pz_summary, pz_summary_file, sep="\t", quote=FALSE, na="NA")
}

if (nrow(qc$se_n_summary) > 0) {
    fwrite(qc$se_n_summary, sen_file, sep="\t", quote=FALSE, na="NA")
}

# ============================================================
# CONVERT CLEANED DATA TO GWAMA FORMAT
#
# GWAMA quantitative input: MARKERNAME, EA, NEA, BETA, SE
# GWAMA binary input:       MARKERNAME, EA, NEA, OR, OR_95L, OR_95U
# N and EAF are included because GWAMA supports them and they are
# useful for meta-analysis summaries. CHR and POS are retained as
# additional marker-location columns.
#
# MARKERNAME uses the standardized key created during QC:
#   CHR:POS:EA:NEA
# ============================================================

if (nrow(clean) > 0) {

    if (trait_type == "quantitative") {

        gwama <- clean[, .(
            MARKERNAME = VARIANT_ID,
            CHR = CHR,
            POS = POS,
            EA = EA,
            NEA = NEA,
            EAF = EAF,
            N = N,
            BETA = BETA,
            SE = SE
        )]

    } else {

        # For binary-trait GWAS from REGENIE/SAIGE, BETA is assumed
        # to be the log odds ratio for the effect allele.
        gwama <- clean[, .(
            MARKERNAME = VARIANT_ID,
            CHR = CHR,
            POS = POS,
            EA = EA,
            NEA = NEA,
            EAF = EAF,
            N = N,
            OR = exp(BETA),
            OR_95L = exp(BETA - 1.96 * SE),
            OR_95U = exp(BETA + 1.96 * SE)
        )]
    }

    # Final GWAMA safety checks.
    if (anyDuplicated(gwama$MARKERNAME)) {
        stop("Duplicated MARKERNAME values remain in GWAMA output.")
    }

    if (any(is.na(gwama$MARKERNAME) | gwama$MARKERNAME == "")) {
        stop("Missing MARKERNAME detected in GWAMA output.")
    }

    if (any(is.na(gwama$EA) | is.na(gwama$NEA))) {
        stop("Missing allele values detected in GWAMA output.")
    }

    if (any(gwama$EAF < 0 | gwama$EAF > 1, na.rm=TRUE)) {
        stop("EAF outside [0,1] detected in GWAMA output.")
    }

    if (trait_type == "quantitative") {
        if (any(!is.finite(gwama$BETA) | !is.finite(gwama$SE) | gwama$SE <= 0)) {
            stop("Invalid BETA/SE detected in quantitative GWAMA output.")
        }
    } else {
        if (any(!is.finite(gwama$OR) | !is.finite(gwama$OR_95L) |
                !is.finite(gwama$OR_95U) | gwama$OR <= 0 |
                gwama$OR_95L <= 0 | gwama$OR_95U <= 0)) {
            stop("Invalid OR or confidence interval detected in binary GWAMA output.")
        }
    }

    gwama_file <- file.path(
        opt$output_dir,
        paste0(output_prefix, ".GWAMA.txt")
    )

    fwrite(
        gwama,
        gwama_file,
        sep="\t",
        quote=FALSE,
        na="NA"
    )

} else {
    gwama <- data.table()
    gwama_file <- NA_character_
}

# ============================================================
# OPTIONAL QC PLOTS
# ============================================================

if (nrow(clean) > 0) {

    png(
        file.path(opt$output_dir, paste0(output_prefix, ".PZ.png")),
        width=1800, height=1800, res=200
    )
    plot_pz(clean, opt$cohort)
    dev.off()

    png(
        file.path(opt$output_dir, paste0(output_prefix, ".QQ.png")),
        width=1800, height=1800, res=200
    )
    plot_qq(clean, opt$cohort)
    dev.off()
}

# ------------------------------------------------------------
# INSERT YOUR EXISTING POST-QC PIPELINE CODE HERE
#
# Use:
#   clean
# or:
#   qc$clean_data
#
# GWAMA-formatted data are available as:
#   gwama
# ------------------------------------------------------------

cat("\n========================================\n")
cat("GWAS QC COMPLETE\n")
cat("========================================\n")
cat("Input variants :", format(nrow(dat), big.mark=","), "\n")
cat("Clean variants :", format(nrow(clean), big.mark=","), "\n")
cat("Clean file     :", clean_file, "\n")
cat("QC summary     :", summary_file, "\n")
cat("QC log         :", log_file, "\n")
cat("GWAMA file     :", ifelse(is.na(gwama_file), "not created", gwama_file), "\n")
cat("========================================\n")
