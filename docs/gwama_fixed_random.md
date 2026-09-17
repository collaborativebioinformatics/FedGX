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
- GWAMA available as `GWAMA` on `PATH`.
- Python 3.10 or newer for the comparison and plotting scripts.
- Matplotlib for Manhattan plots.

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
python3 -m pip install matplotlib
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

`N` and `EAF` are recommended when available. For Manhattan plots,
`MARKERNAME` must begin with chromosome and position:

```text
1:123456:A:C
```

Files created by `Scripts/gwas_cohort_qc_with_gwama.R` retain this coordinate
information. Identifiers such as `rs123` or `SNP1` can be meta-analysed but
cannot be positioned on a Manhattan plot without a separate map.

The runner does not require a particular input directory. Every site file is
passed as an explicit command-line path. `Scripts/gwas_cohort_qc_with_gwama.R`
writes `<output-prefix>.GWAMA.txt.gz` into `--output-dir`; when `--output-dir`
is omitted, it writes to the current directory.

A recommended central layout is:

```text
runs/<phenotype>/<method>/
  inputs/
    site1.GWAMA.txt.gz
    site2.GWAMA.txt.gz
    site3.GWAMA.txt.gz
  meta.fixed.out
  meta.random.out
  meta.comparison.tsv
```

For example, use `--output-dir runs/phenotype/regenie/inputs` and
`--output-prefix site1` when creating the first site's GWAMA file, or copy the
site-produced summary file to that location on the central server.

## Run fixed and random models

Run commands from the repository root. For a binary/case-control trait:

```bash
mkdir -p runs/phenotype/regenie/inputs

bash Scripts/run_gwama.sh or --model both \
  runs/phenotype/regenie/meta \
  runs/phenotype/regenie/inputs/site1.GWAMA.txt.gz \
  runs/phenotype/regenie/inputs/site2.GWAMA.txt.gz \
  runs/phenotype/regenie/inputs/site3.GWAMA.txt.gz
```

For a quantitative trait, replace `or` with `qt`:

```bash
bash Scripts/run_gwama.sh qt --model both \
  runs/phenotype/regenie/meta \
  runs/phenotype/regenie/inputs/site1.GWAMA.txt.gz \
  runs/phenotype/regenie/inputs/site2.GWAMA.txt.gz \
  runs/phenotype/regenie/inputs/site3.GWAMA.txt.gz
```

To run only one model, use `--model fixed` or `--model random`. Omitting
`--model` preserves the original fixed-effect behavior:

```bash
bash Scripts/run_gwama.sh or runs/phenotype/regenie/meta \
  runs/phenotype/regenie/inputs/site1.GWAMA.txt.gz \
  runs/phenotype/regenie/inputs/site2.GWAMA.txt.gz
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

### `GWAMA executable not found in PATH`

Add the directory containing the executable, not the executable itself:

```bash
export PATH="/path/to/GWAMA-directory:$PATH"
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
