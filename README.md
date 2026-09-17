# FedGX: Federated GWAS Across Cohort Sites

**Federated learning software for multi-tool genome-wide association studies (GWAS) across distributed cohort sites.**
Built at the Nordic Biobank × NVIDIA Hackathon.

![Logo](docs/logo.jpeg)

![Status](https://img.shields.io/badge/status-hackathon_prototype-orange?style=for-the-badge)
![License](https://img.shields.io/badge/license-MIT-blue?style=for-the-badge)
![FL Framework](https://img.shields.io/badge/NVIDIA%20FLARE-2.7.1-76B900?style=for-the-badge&logo=nvidia&logoColor=white)
![Sites](https://img.shields.io/badge/sites-3-9146FF?style=for-the-badge)
![Python](https://img.shields.io/badge/python-3.10%2B-3776AB?style=for-the-badge&logo=python&logoColor=white)
![R](https://img.shields.io/badge/R-GWAMA%20QC-276DC3?style=for-the-badge&logo=r&logoColor=white)

---

## Table of Contents

- [Summary](#summary)
- [Architecture](#architecture)
- [Quickstart](#quickstart)
- [Synthetic Data Specifications](#synthetic-data-specifications)
- [Detailed Setup Instructions](#detailed-setup-instructions)
- [Technologies](#technologies)
- [References](#references)
- [License](#license)
- [Team](#team)

---

## Summary

At the Nordic Biobank × NVIDIA Hackathon, we built software for running a federated GWAS across three cohort sites **without centralizing raw genotype data**. Each site runs its own association analysis locally; only summary statistics are shared and combined centrally via meta-analysis.

Development and testing were carried out on [HUNT Cloud](https://www.ntnu.edu/mh/hunt/hunt-cloud) and across two [Brev](https://brev.dev) GPU instances.

This project extends [FedGen](https://github.com/collaborativebioinformatics/FedGen), building on its federated-learning scaffolding with a synthetic multi-ancestry data pipeline, a REGENIE-based per-site GWAS workflow, and a fixed-/random-effects GWAMA meta-analysis with automated comparison and plotting.

> [!TIP]
> No real cohort data on hand? Skip straight to [Synthetic Data Specifications](#synthetic-data-specifications) — `Scripts/generateData.sh` builds a full 3-site multi-ancestry dataset locally.

---

## Architecture

![Federated GWAS architecture](docs/FederationFigure_MR.png)

---

## Quickstart

### 1. Start the NVFLARE Dashboard and FL Server

*(Server-side setup — see your NVFLARE dashboard documentation.)*

### 2. Start an NVFLARE Client on Brev / HUNT Cloud

#### 2.1 Create a GPU instance on Brev

On the **Brev website**:

- Create **one GPU instance per site**
- Example configuration:

  | Setting | Value |
  |---|---|
  | Name | `site1` |
  | GPU | 1× NVIDIA L4 |
  | CPU | 16 cores |
  | RAM | 64 GB |

#### 2.2 Connect to the instance

```bash
brev shell site1
```

Use a terminal multiplexer so the session survives disconnects (optional but recommended):

```bash
tmux new -s nvflare
```

#### 2.3 Set up the Python environment

```bash
python3 -m venv venv_nvflare
source venv_nvflare/bin/activate

pip install nvflare[PT] torch torchvision tensorboard
```

Verify the install:

```bash
nvflare --version
```

### 3. Copy and start the NVFLARE client startup kit

#### 3.1 Copy the client kit from your local machine

On your **local machine**:

```bash
brev copy <local_path_to_client_kit> site1:<remote_path>
```

On the **Brev instance**:

```bash
sudo apt update
sudo apt install -y unzip

unzip -d <client_name> -P <PIN> <client_kit.zip>
cd <client_name>
```

#### 3.2 Start the NVFLARE client

```bash
./startup/start.sh
```

Check the logs to confirm the client connected to the NVFLARE server/dashboard.

### 4. Clone the FedGX repository

```bash
git clone https://github.com/collaborativebioinformatics/FedGX
chmod +x FedGen/scripts/*.sh
```

### 5. Get site data

Two options, depending on your setup:

- **Real cohort data:** download your assigned site from S3 (see your site coordinator for credentials and the bucket path).
- **Synthetic data:** generate a full synthetic multi-ancestry dataset locally with `Scripts/generateData.sh` — no S3 access needed. See [Synthetic Data Specifications](#synthetic-data-specifications) below for exactly what it produces.

> [!WARNING]
> The S3 download command/script path for real cohort data needs confirming — the previous reference to `scripts/download_site_from_s3.sh` no longer matches a file in this repo.

### 6. Run REGENIE per site (outside NVFLARE)

Run REGENIE independently per site (not through NVFLARE) to verify all dependencies are working before attempting a federated run:

```bash
cd ~/data
./../FedGen/scripts/run_regenie_site.sh <siteNumber>
```

Monitor the logs and outputs to confirm successful completion.

**Runtime:** ~30–45 minutes total

- Step 1 (LOCO model): 15–30 min
- Step 2 (association testing): 10–20 min

### 7. Run the federated GWAS job (NVFLARE)

Instead of running REGENIE independently on each site and manually aggregating results, submit a federated GWAS job that automates the whole workflow across all sites using NVIDIA FLARE.

The federated job handles:

- Distributing analysis scripts to all clients
- Running local GWAS analysis with REGENIE on each site
- Collecting summary statistics from all sites
- Performing meta-analysis with GWAMA on the server

> [!WARNING]
> `jobs/fed_gwas/README.md` is referenced here in the previous draft but doesn't exist in this repo yet — add it, or point this section at wherever the federated job config actually lives.

### 8. Run the central GWAMA meta-analysis

The central runner supports fixed-effect (`fixed`), random-effects (`random`), or paired (`both`) GWAMA analyses for binary and quantitative traits. Paired runs also produce a marker-level FFX/RFX comparison, and a separate script creates matched fixed- and random-effects Manhattan plots.

Before the central run, each site uses `Scripts/gwas_cohort_qc_with_gwama.R` to QC its REGENIE, SAIGE, or PLINK summary results and produce `<site>.GWAMA.txt.gz`. This requires an explicit mapping from the source columns to chromosome, position, alleles, frequency, effect, standard error, sample size, and P-value — see [`docs/gwama_fixed_random.md`](docs/gwama_fixed_random.md) for a complete REGENIE example and the requirements for other formats.

Quick start for a binary trait:

```bash
mkdir -p runs/phenotype/regenie/inputs

bash Scripts/run_gwama.sh or --model both \
  runs/phenotype/regenie/meta \
  runs/phenotype/regenie/inputs/site1.GWAMA.txt.gz \
  runs/phenotype/regenie/inputs/site2.GWAMA.txt.gz \
  runs/phenotype/regenie/inputs/site3.GWAMA.txt.gz

python3 Scripts/plot_gwama_manhattan.py \
  --fixed runs/phenotype/regenie/meta.fixed.out \
  --random runs/phenotype/regenie/meta.random.out \
  --output-prefix runs/phenotype/regenie/meta
```

> [!IMPORTANT]
> Run REGENIE, SAIGE, and PLINK results as **separate** meta-analyses — don't mix methods in one GWAMA run. All contributing sites must share the same phenotype, trait definition, and genome build.

See [`docs/gwama_fixed_random.md`](docs/gwama_fixed_random.md) for installation checks, input requirements, the complete command reference, output formats, result interpretation, testing, and troubleshooting.

### 9. Notes & best practices

- Use **one Brev instance per NVFLARE client**
- Always run the NVFLARE client inside a virtual environment
- Prefer **IAM roles** over static AWS credentials
- Validate GPU availability: `nvidia-smi`
- Use `tmux` or `screen` to keep long-running jobs alive

---

## Synthetic Data Specifications

Generated by `Scripts/generateData.sh` — a full run needs no real cohort data or S3 access.

![Samples](https://img.shields.io/badge/samples-88K--110K%2Fsite-6A4C93?style=flat-square)
![SNPs](https://img.shields.io/badge/SNPs-450K--520K%2Fsite-1982C4?style=flat-square)
![Heritability](https://img.shields.io/badge/h²-0.20--0.25-8AC926?style=flat-square)
![Ancestries](https://img.shields.io/badge/ancestries-EUR%20%7C%20EAS%20%7C%20AFR-FFCA3A?style=flat-square)
![Causals](https://img.shields.io/badge/causal%20variants-20%2Fsite-FF595E?style=flat-square)

### Genotypes

| Parameter | Value | Notes |
|---|---|---|
| Format | PLINK binary (`.bed`/`.bim`/`.fam`) | Standard genetic format |
| Variants per site | 450K–520K SNPs | Master grid: 520K; sites subset to their own count |
| Samples per site | 88K–110K individuals | Site 1: 100K, Site 2: 95K, Site 3: 110K |
| Chromosomes | 22 autosomes | hg38 assumed (implicit in LDAK) |
| Build | hg38 | Not explicit in the script; inferred from the PLINK standard |
| MAF range | 0.01–0.50 (ancestry-dependent) | EUR: 0.01–0.50; EAS: 0.01–0.45; AFR: 0.005–0.50 |
| LD structure | Realistic | Generated by LDAK; reflects linkage disequilibrium |
| Population structure | 1–3 subpopulations per site | Site 1: 1 (homogeneous); Site 2: 2; Site 3: 3 |

### Phenotype

| Parameter | Value | Notes |
|---|---|---|
| File format | Space-delimited (`FID IID Pheno`) | Standard PLINK phenotype format |
| Trait | Parkinson's disease (binary) | 0 = control, 1 = case |
| Prevalence | 1% | Realistic for elderly populations |
| Heritability (h²) | 0.20–0.25 on liability scale | Site 1: 0.25; Site 2: 0.22; Site 3: 0.20 |
| Causal variants | 20 per site | 12 shared (identical across sites) + 8 site-specific |
| Effect size model | LDAK-Thin (power = −0.25) | Realistic allelic architecture |
| Effect scale | 1.00, 0.85, 0.70 across sites | Applied to the first half of shared causals; produces heterogeneity |
| Fixed effects | Yes (`USE_FIXED_EFFECTS=1`) | Deterministic, reproducible effects; first half scaled by ancestry |

### Covariates

| Parameter | Value | Notes |
|---|---|---|
| File | `site{N}_geno.covar` | Auto-generated by LDAK `--make-snps` |
| Variables | Age, sex, principal components | Exact covariates depend on the LDAK output structure |
| Variance explained | ~10% of phenotypic variation (`COVAR_HER=0.1`) | Controls confounding; realistic level |

### Shared causal architecture

| Parameter | Value | Notes |
|---|---|---|
| Shared causal variants | 12 | Identical positions across all three sites |
| Selection method | Fixed fractional positions in the genome | Every site generates 520K variants; same physical locations |
| Effect sizes (shared) | ±(0.60–0.84) | Deterministic function of position; alternating sign |
| Ancestry scaling | 1.00 (EUR), 0.85 (EAS), 0.70 (AFR) | Applied only to the first 6 shared causals to preserve detectability |

### Site-specific causal architecture

| Parameter | Value | Notes |
|---|---|---|
| Site-specific causals per site | 8 | Disjoint between sites; non-overlapping variants |
| Selection method | Deterministic band sampling (offset by site) | Spreads variants across the genome; ensures no overlap |
| Effect sizes (site-specific) | ±(0.70–0.91) | Larger than shared causals; produces local signal |
| Total causals per site | 20 | 12 shared + 8 specific |

### Key design decisions

1. **Master grid approach** — all sites generate the same 520K-variant skeleton, then subset to their own count. This puts shared causal variants at identical physical positions without any coordination step.
2. **Ancestry heterogeneity** — MAF and effect scales differ per site, producing realistic I² (heterogeneity) values in meta-analysis while keeping shared signals detectable.
3. **Fixed effects** — deterministic and reproducible across runs, which makes validation and debugging tractable.
4. **Power parameter** — −0.25 reflects realistic allelic architecture, where common variants have smaller effects than rare variants.

---

## Detailed Setup Instructions

### Prerequisites

> [!NOTE]
> List required tools and versions here (Python, R, PLINK, REGENIE, GWAMA, LDAK, Docker, AWS CLI, NVFLARE) with install commands or links.

### Download workflow

```bash
# 1. Clone the repository (if not already done)
git clone https://github.com/collaborativebioinformatics/FedGen.git
cd FedGen

# 2. Download your assigned site (e.g., Site 3)
./scripts/download_site_from_s3.sh 3

# 3. Verify the download
ls -lh data/simulated_sites/site3/
# Should show ~15 GB total:
# - site3_geno.bed   (~12-13 GB)
# - site3_geno.bim   (~10-20 MB)
# - site3_geno.fam   (~2-3 MB)
# - site3_pheno.pheno (~2 MB)
# - site3_geno.covar (~5-10 MB)
```

> [!WARNING]
> This points at the `FedGen` repo's `download_site_from_s3.sh`, which doesn't exist in `FedGX`/`Scripts/`. Confirm whether real-data users should be pointed at FedGen instead, or whether this script needs to be ported into this repo.

### REGENIE analysis workflow

> [!NOTE]
> **@Allan:** describe the REGENIE step-1/step-2 commands used per site, required flags, and expected runtime here.

#### Run analysis

> [!NOTE]
> Command block for running the analysis.

### Understanding results

> [!NOTE]
> **@Xiaoping:** explain how to interpret REGENIE/GWAMA output for this pipeline.

#### Output files

> [!NOTE]
> List the output files and what each contains.

#### Association results format

> [!NOTE]
> Document the column layout of the association output file.

**Key columns:**

> [!NOTE]
> List and describe the key columns (e.g., CHR, POS, A1, A1_FREQ, BETA, SE, P).

#### Finding genome-wide significant hits

> [!NOTE]
> Describe the filtering convention (e.g., `P < 5e-8`) and how to extract hits.

#### Manhattan plot (R)

> [!NOTE]
> R plotting instructions, or a pointer to `Scripts/plot_gwama_manhattan.py` if the R version was dropped in favor of the Python one.

---

## Technologies

| Tool | Role |
|---|---|
| [LDAK](https://dougspeed.com/) v6.1 | Data generation |
| [REGENIE](https://rgcgithub.github.io/regenie/) v4.1 | GWAS analysis |
| [GWAMA](https://genomics.ut.ee/en/tools/gwama) | Meta-analysis |
| Docker | Containerization |
| AWS S3 | Data storage |
| NVIDIA FLARE 2.7.1 | Federated learning framework |
| Brev | Compute for distributed sites |

---

## References

### Software citations

- **LDAK:** Speed et al. (2020). Improved heritability estimation from genome-wide SNPs. *Nature Genetics*. https://doi.org/10.1038/s41588-019-0530-8
- **REGENIE:** Mbatchou et al. (2021). Computationally efficient whole-genome regression for quantitative and binary traits. *Nature Genetics*. https://doi.org/10.1038/s41588-021-00870-7
- **GWAMA:** Mägi et al. (2010). GWAMA: software for genome-wide association meta-analysis. *BMC Bioinformatics*. https://doi.org/10.1186/1471-2105-11-288
- **This project:** *(add citation when published)*

### Documentation links

- [NVFLARE Documentation](https://nvflare.readthedocs.io/)
- [FedGen Repository](https://github.com/collaborativebioinformatics/FedGen)
- [Brev Platform](https://brev.dev)
- [REGENIE Documentation](https://rgcgithub.github.io/regenie/)
- [LDAK Documentation](https://dougspeed.com/)
- [PLINK File Formats](https://www.cog-genomics.org/plink/1.9/formats)

---

## License

Data and scripts: MIT License (see repository root).

---

## Team

| Contributor | Sections | Affiliation |
|---|---|---|
| Marlene Rietz | 1–3 | Steno Diabetes Center Odense, Odense, Denmark; P1 Pioneer Center for Artificial Intelligence, University of Copenhagen, Copenhagen, Denmark; Department of Laboratory Medicine, Karolinska Institutet, Stockholm, Sweden |
| Allan Lind-Thomsen | 4 | *(TODO: affiliation)* |
| Xiaoping Wu | 5 | Department of Obstetrics and Gynecology, Institute of Clinical Sciences, Sahlgrenska Academy, University of Gothenburg, Gothenburg, Sweden |
| Moh Sallam | 6 | Center for Quantitative Genetics and Genomics and Pioneer Center for Smartbiomed, Aarhus University |
| Pravesh Parekh | 7–8 | J. Craig Venter Institute, San Diego, California, USA; Centre for Precision Psychiatry, University of Oslo, Oslo, Norway |

---

<p align="center"><sub>Built at the Nordic Biobank × NVIDIA Hackathon</sub></p>
