# Central GWAMA fixed- and random-effects analysis

This guide describes the central meta-analysis stage after each participating
site has completed GWAS quality control and converted its summary statistics to
GWAMA format.

The workflow runs each GWAS method separately:

```text
QCed site results
  -> fixed-effect GWAMA
  -> random-effects GWAMA
  -> FFX/RFX comparison table
  -> fixed and random Manhattan plots
```

Do not combine REGENIE, SAIGE, and PLINK site results in one analysis. Use a
separate output directory for each method.

## Requirements

- Linux or WSL with Bash.
- GWAMA available as `GWAMA` on `PATH`, at
  `$FEDGX_TOOLS_ROOT/GWAMA`, or through the exact `$FEDGX_GWAMA_BIN` path.
- R with the `optparse` and `data.table` packages for site-level QC and
  conversion. R is not required if GWAMA-formatted files already exist.
- Python 3.10 or newer for the comparison and plotting scripts.
- Matplotlib for Manhattan plots and pandas for the optional dashboard builder.

Check the environment:

```bash
command -v GWAMA
GWAMA --help | head
python3 --version
python3 -c "import matplotlib; print(matplotlib.__version__)"
```

If Matplotlib is missing, activate the project virtual environment and install
it there:

```bash
python3 -m pip install matplotlib pandas
```

## Input contract

Provide at least two site-level GWAMA files from the same analysis. Before
combining them, confirm that all sites use:

- The same phenotype and phenotype coding.
- The same binary or quantitative trait type.
- The same genome build, for example GRCh38.
- The same GWAS method and broadly equivalent covariates.
- Harmonized effect and non-effect alleles.

Binary-trait files require these core columns:

```text
MARKERNAME EA NEA OR OR_95L OR_95U
```

Quantitative-trait files require:

```text
MARKERNAME EA NEA BETA SE
```

`N` and `EAF` are required by the current runner. For Manhattan plots,
`MARKERNAME` must begin with chromosome and position:

```text
1:123456:A:C
```

Files created by `Scripts/gwas_cohort_qc_with_gwama.R` retain this coordinate
information. Identifiers such as `rs123` or `SNP1` can be meta-analysed but
cannot be positioned on a Manhattan plot without a separate map.

The runner does not require a particular input directory. Every site file is
passed as an explicit command-line path. `Scripts/gwas_cohort_qc_with_gwama.R`
writes `<output-prefix>.GWAMA.txt` into `--output-dir`; when `--output-dir`
is omitted, it writes to the current directory.

A recommended central layout is:

```text
runs/<phenotype>/<method>/
  inputs/
    site1.GWAMA.txt
    site2.GWAMA.txt
    site3.GWAMA.txt
  meta.fixed.out
  meta.random.out
  meta.comparison.tsv
```

For example, use `--output-dir runs/phenotype/regenie/inputs` and
`--output-prefix site1` when creating the first site's GWAMA file, or copy the
site-produced summary file to that location on the central server.

## Create a GWAMA file at each site

Run `Scripts/gwas_cohort_qc_with_gwama.R` at each participating site before
central meta-analysis. It performs cohort-level QC and writes an uncompressed
GWAMA input file without transferring individual-level genotype or phenotype
data.

The converter does not guess the source software. Inspect the GWAS header and
map each required concept to its actual column name:

| Concept | Command option | REGENIE example |
| --- | --- | --- |
| Chromosome | `--col-chr` | `CHROM` |
| Position | `--col-pos` | `GENPOS` |
| Variant ID | `--col-id` | `ID` |
| Effect allele | `--col-ea` | `ALLELE1` |
| Non-effect allele | `--col-nea` | `ALLELE0` |
| Effect-allele frequency | `--col-eaf` | `A1FREQ` |
| Effect or log odds ratio | `--col-beta` | `BETA` |
| Standard error | `--col-se` | `SE` |
| Per-variant sample size | `--col-n` | `N` |
| P-value or -log10(P) | `--col-p` or `--col-log10p` | `LOG10P` |

For a quantitative REGENIE result, an example site command is:

```bash
mkdir -p site_outputs/phenotype/regenie

Rscript Scripts/gwas_cohort_qc_with_gwama.R \
  --input results/regenie_step2_Phen1.regenie.gz \
  --cohort SITE1 \
  --build GRCh38 \
  --trait-type quantitative \
  --info-threshold 0.30 \
  --min-n 30 \
  --min-mac 6 \
  --remove-palindromic TRUE \
  --autosomes-only TRUE \
  --snp-only TRUE \
  --remove-duplicates TRUE \
  --filter-test TRUE \
  --test-value ADD \
  --col-chr CHROM \
  --col-pos GENPOS \
  --col-id ID \
  --col-ea ALLELE1 \
  --col-nea ALLELE0 \
  --col-eaf A1FREQ \
  --col-beta BETA \
  --col-se SE \
  --col-n N \
  --col-info INFO \
  --col-log10p LOG10P \
  --col-test TEST \
  --col-chisq CHISQ \
  --output-dir site_outputs/phenotype/regenie \
  --output-prefix site1
```

The central input produced by this example is:

```text
site_outputs/phenotype/regenie/site1.GWAMA.txt
```

Repeat with `site2`, `site3`, and so on. Transfer only the approved summary
outputs to the central server, then place them in the chosen central input
directory.

For a binary trait, set `--trait-type binary`. The column passed to
`--col-beta` must contain a log odds ratio; the converter calculates OR and its
95% confidence interval. Do not pass an untransformed OR as `--col-beta`.

For SAIGE and PLINK, use the same command but replace the column mappings with
the actual headers in those outputs. The input must include effect size (BETA,
or log OR for a binary trait), SE, EAF, and per-variant N. If a required value
is absent, add it during an upstream standardization step rather than assigning
an unrelated column.

Besides the `.GWAMA.txt` file, the site script writes cleaned statistics,
QC summaries and logs, P-Z and SE-N summaries, and cohort-level P-Z and QQ
plots. Review the QC log before releasing the GWAMA file for central analysis.

## Run fixed and random models

Run commands from the repository root. For a binary/case-control trait:

```bash
mkdir -p runs/phenotype/regenie/inputs

bash Scripts/run_gwama.sh or --model both \
  runs/phenotype/regenie/meta \
  runs/phenotype/regenie/inputs/site1.GWAMA.txt \
  runs/phenotype/regenie/inputs/site2.GWAMA.txt \
  runs/phenotype/regenie/inputs/site3.GWAMA.txt
```

For a quantitative trait, replace `or` with `qt`:

```bash
bash Scripts/run_gwama.sh qt --model both \
  runs/phenotype/regenie/meta \
  runs/phenotype/regenie/inputs/site1.GWAMA.txt \
  runs/phenotype/regenie/inputs/site2.GWAMA.txt \
  runs/phenotype/regenie/inputs/site3.GWAMA.txt
```

To run only one model, use `--model fixed` or `--model random`. Omitting
`--model` preserves the original fixed-effect behavior:

```bash
bash Scripts/run_gwama.sh or runs/phenotype/regenie/meta \
  runs/phenotype/regenie/inputs/site1.GWAMA.txt \
  runs/phenotype/regenie/inputs/site2.GWAMA.txt
```

## Outputs

With `--model both`, the output prefix above produces:

| File | Description |
| --- | --- |
| `meta.gwama.in` | Deterministic list of contributing site files. |
| `meta.fixed.out` | GWAMA fixed-effect association results. |
| `meta.fixed.log.out` | Fixed-effect GWAMA log. |
| `meta.fixed.err.out` | Fixed-effect GWAMA error stream. |
| `meta.random.out` | GWAMA random-effects association results. |
| `meta.random.log.out` | Random-effects GWAMA log. |
| `meta.random.err.out` | Random-effects GWAMA error stream. |
| `meta.comparison.tsv` | Marker-level FFX/RFX comparison. |

The comparison table contains fixed and random effects, standard errors and
P-values, absolute differences, Cochran's Q, the heterogeneity P-value, I2,
the number of studies, and marker-matching status.

When I2 is zero, the random-effects result may be identical to the fixed-effect
result. With substantial heterogeneity, RFX commonly has a larger standard
error and a less significant P-value.

## Run through the federated server controller

`Scripts/serverSide_job.py` is the integrated NVFLARE entry point. It packages
the site driver and QC converter for clients and packages the importable
aggregator, GWAMA runner, comparison helper, and selected post-processing
helpers for the server. Run it from any directory; packaged source paths are
resolved relative to the script itself.

```bash
python3 Scripts/serverSide_job.py \
  --env prod \
  --n_clients 3 \
  --method regenie \
  --trait_type binary \
  --model both \
  --tools_root /home/ubuntu/tools \
  --startup_kit /path/to/server/startup-kit \
  --username your-nvflare-user
```

The default `--model both --plots` combination runs fixed and random GWAMA,
checks that their expected files exist, creates the comparison table, and
creates both Manhattan PNGs. Use `--dashboard` to add separate interactive
fixed- and random-effects dashboard directories. Use `--no-plots` only when
deliberately running a single model.

## Create the Manhattan plots

Generate two PNG files from a paired run:

```bash
python3 Scripts/plot_gwama_manhattan.py \
  --fixed runs/phenotype/regenie/meta.fixed.out \
  --random runs/phenotype/regenie/meta.random.out \
  --output-prefix runs/phenotype/regenie/meta
```

Outputs:

```text
runs/phenotype/regenie/meta.fixed.manhattan.png
runs/phenotype/regenie/meta.random.manhattan.png
```

The two plots use the same chromosome layout and Y-axis scale, allowing direct
visual comparison. Both include the default genome-wide significance line at
`P = 5e-8`.

Change the significance threshold or image resolution when needed:

```bash
python3 Scripts/plot_gwama_manhattan.py \
  --fixed runs/phenotype/regenie/meta.fixed.out \
  --random runs/phenotype/regenie/meta.random.out \
  --output-prefix runs/phenotype/regenie/meta \
  --threshold 5e-8 \
  --dpi 300
```

## Build the optional interactive dashboard

The dashboard is an additional reporting layer, not a replacement for the two
static Manhattan PNGs. It accepts the same central `site*.GWAMA.txt` files and
one fixed- or random-effects GWAMA output. Binary OR/confidence-interval values
are converted back to log effects for the forest plot.

```bash
python3 Scripts/build_fedx_dashboard.py \
  --site runs/phenotype/regenie/inputs/site1.GWAMA.txt \
  --site runs/phenotype/regenie/inputs/site2.GWAMA.txt \
  --site runs/phenotype/regenie/inputs/site3.GWAMA.txt \
  --meta runs/phenotype/regenie/meta.fixed.out \
  --model fixed \
  --trait-type binary \
  --method regenie \
  --html-template Scripts/fedx_dashboard.html \
  --out-dir runs/phenotype/regenie/dashboard/fixed
```

Open `runs/phenotype/regenie/dashboard/fixed/index.html`. The generated
`fedx_data.js` contains site-level summary statistics, so treat the dashboard
directory as analysis output and do not publish it without the same approval
used for the underlying summary statistics. To keep a full-GWAS browser view
responsive, the builder embeds at most 50,000 variants by default, retaining
the strongest associations and genome-wide coverage. Use `--max-variants 0`
only when you intentionally want every variant; full-data QQ plots from the R
QC step remain the authoritative diagnostics when the dashboard is sampled.

## Validate an installation

Run the automated comparison and plotting tests:

```bash
python3 -B -m unittest discover -s Scripts/tests -v
bash -n Scripts/run_gwama.sh
```

The test suite does not replace an end-to-end run with the installed GWAMA
binary. Before analysing real data, run a small multi-site example through
`--model both` and confirm that both `.err.out` files are empty.

## Troubleshooting

### `GWAMA was not found`

Choose one supported resolution method:

```bash
export PATH="/path/to/GWAMA-directory:$PATH"
# or
export FEDGX_TOOLS_ROOT="/path/to/tools"
# or
export FEDGX_GWAMA_BIN="/exact/path/to/GWAMA"
```

### Python reports version 2.7

Activate a Python 3 virtual environment and clear Bash's command cache:

```bash
source .venv/bin/activate
hash -r
python3 --version
```

### `syntax error near unexpected token $'{\r''`

The shell script has Windows CRLF line endings. A normal Git checkout respects
the repository `.gitattributes` file. For an individually copied file, convert
it on Linux:

```bash
sed -i 's/\r$//' Scripts/run_gwama.sh
```

### Fixed and random results are identical

Check `q_p-value` and `i2`. Identical estimates are expected when GWAMA detects
little or no heterogeneity.

### `no plottable variants`

Manhattan plots require marker names beginning with `CHR:POS`. Recreate the
GWAMA input with coordinate-based marker identifiers or provide a future map
step before plotting.

## Reproducibility notes

- Keep generated `runs/`, GWAMA `.out`, log, and PNG files out of Git.
- Record the input file names, genome build, phenotype, GWAS method, and GWAMA
  version for every central run.
- Preserve site order and archive `meta.gwama.in` with the analysis results.
- Review allele, effect, and strand warnings before interpreting associations.
