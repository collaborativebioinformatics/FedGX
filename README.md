# FedGX: Federated GWAS Across Cohort Sites

**Federated learning software for multi-tool genome-wide association studies (GWAS) across distributed cohort sites.** Built at the Nordic Biobank × NVIDIA Hackathon.

![Logo](docs/logo.jpeg)

---

## The problem

- Genome-wide association studies need **large, diverse cohorts** to find real signal
- But cohort data usually can't leave the site it came from — privacy, governance, legal restrictions
- Result: most GWAS run on a single biobank, missing ancestries it doesn't have

---

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

Development and testing were carried out on [HUNT Cloud](https://www.ntnu.edu/mh/hunt/hunt-cloud) and across two
[Brev](https://brev.dev) GPU instances.

This project extends
[FedGen](https://github.com/collaborativebioinformatics/FedGen), building on its federated-learning scaffolding with a synthetic multi-ancestry data pipeline, a REGENIE-based per-site GWAS workflow, and a fixed-/random-effects GWAMA meta-analysis with automated comparison and plotting.

> [!TIP] No real cohort data on hand? Skip straight to [Synthetic Data
> Specifications](#synthetic-data-specifications) —
> `Scripts/generateData.sh` builds a full 3-site multi-ancestry dataset
> locally.

---

## Architecture
![Federated GWAS architecture](docs/FederationFigure_MR.png)

---

## Quickstart
### 1. Start the NVFLARE Dashboard and FL Server

To find information on how to set up the NVFLARE dashboard, see this link: <https://nvflare.readthedocs.io/en/2.4/real_world_fl/workspace.html>

### 2. Start all client instances, e.g., secure servers hosting cohort data

Connect to each site's instance (Brev, HUNT Cloud, or your own infrastructure — the connection method depends on where that site is hosted) and set up the Python environment:

``` bash
python3 -m venv venv_nvflare
source venv_nvflare/bin/activate

pip install nvflare[PT] torch torchvision tensorboard pandas matplotlib
```

Verify the install:

``` bash
nvflare --version
```

### 3. Copy and start the NVFLARE client startup kit

The client startup kit is a PIN-protected zip you download once per client from the NVFLARE dashboard/admin console. Besides the start.sh script, it bundles the set of SSL certificates that let this client authenticate to the NVFLARE server: a client certificate and private key, plus the server's root CA certificate, so the client and server can establish mutual TLS. You only need to redo this step the first time you set up a client, or if the kit is later regenerated (e.g. the certificates expire or the client is re-provisioned).

#### 3.1 Copy the client kit from your local machine

On your **local machine** (Brev example — adjust for other hosting):

``` bash
brev copy <local_path_to_client_kit> site1:<remote_path>
```

On the client instance, unzip it with the PIN provided when the kit was created:

``` bash
sudo apt update
sudo apt install -y unzip

unzip -d <client_name> -P <PIN> <client_kit.zip>
cd <client_name>
```

#### 3.2 Start the NVFLARE client

``` bash
./startup/start.sh
```

Check the logs to confirm the client connected to the NVFLARE server/dashboard.

### 4. Clone the FedGX repository

``` bash
git clone https://github.com/collaborativebioinformatics/FedGX
chmod +x FedGX/Scripts/*.sh
```

### 5. Get site data

Two options, depending on your setup:

- **Real cohort data** 
- **Synthetic data** generate a full synthetic multi-ancestry dataset locally with `Scripts/generateData.sh` — no S3 access needed. See
    [Synthetic Data Specifications](#synthetic-data-specifications) below for exactly what it produces.

Configure the local site client to obtain metadata: create `~/fedgx_config/config.json` with the keys `Scripts/local_script_start_gwas.sh` reads at startup:

``` json
{
  "program": "regenie",
  "populationid": "site1",
  "traittype": "binary",
  "build": "GRCh38"
}
```

`traittype` defaults to `binary` and `build` to `GRCh38` if omitted. The REGENIE binary path is set separately via the `REGENIE` environment variable, not through this file.
`programs` should be in the PATH on the client

### 6. Run REGENIE per site (outside NVFLARE)

Run REGENIE independently per site (not through NVFLARE) to verify all dependencies are working before attempting a federated run:

``` bash
cd ~/data
./../FedGX/Scripts/local_script_start_gwas.sh
```

This reads `~/fedgx_config/config.json` (see [Get site data](#5-get-site-data) above), runs REGENIE steps 1–2 against the LDAK-generated genotypes, and writes GWAMA-format output.

Monitor the logs and outputs to confirm successful completion.

**Runtime:** \~30–45 minutes total

- Step 1 (LOCO model): 15–30 min
- Step 2 (association testing): 10–20 min

### 7. Run the federated GWAS job (NVFLARE)

Instead of running REGENIE independently on each site and manually aggregating results, submit a federated GWAS job that automates the whole workflow across all sites using NVIDIA FLARE.

The federated job handles:

- Distributing analysis scripts to all clients
- Running local GWAS analysis with REGENIE on each site
- Collecting summary statistics from all sites
- Performing meta-analysis with GWAMA on the server

The server entry point coordinates the complete workflow. A paired run is the default: fixed effects are the primary result and random effects are produced as a sensitivity analysis. Add `--dashboard` to build optional interactive fixed- and random-effects views after the standard result files and PNGs.

``` bash
python3 Scripts/serverSide_job.py \
  --env prod \
  --n_clients 3 \
  --method regenie \
  --trait_type binary \
  --model both \
  --tools_root /home/ubuntu/tools \
  --startup_kit /path/to/server/startup-kit \
  --username your-nvflare-user \
  --dashboard
```

Without `--dashboard`, the central job still writes the GWAMA fixed/random results, their comparison table, and two matched Manhattan PNGs. The server must have executable `/home/ubuntu/tools/GWAMA`; alternatively set `FEDGX_GWAMA_BIN` to its exact path.

### 8. Run the central GWAMA meta-analysis

The central runner supports fixed-effect (`fixed`), random-effects (`random`), or paired (`both`) GWAMA analyses for binary and quantitative traits. Paired federated runs produce a marker-level FFX/RFX comparison and matched fixed- and random-effects Manhattan plots automatically. The plotting script remains available for standalone reruns.

Before the central run, each site uses `Scripts/gwas_cohort_qc_with_gwama.R` to QC its REGENIE summary results and produce an uncompressed `<site>.GWAMA.txt`. This requires an explicit mapping from the source columns to chromosome, position, alleles, frequency, effect, standard error, sample size, and P-value — see [`docs/gwama_fixed_random.md`](docs/gwama_fixed_random.md) for the complete format requirements. SAIGE and PLINK should only be advertised after their site-driver branches are implemented and tested.

Quick start for a binary trait:

``` bash
mkdir -p runs/phenotype/regenie/inputs

bash Scripts/run_gwama.sh or --model both \
  runs/phenotype/regenie/meta \
  runs/phenotype/regenie/inputs/site1.GWAMA.txt \
  runs/phenotype/regenie/inputs/site2.GWAMA.txt \
  runs/phenotype/regenie/inputs/site3.GWAMA.txt

python3 Scripts/plot_gwama_manhattan.py \
  --fixed runs/phenotype/regenie/meta.fixed.out \
  --random runs/phenotype/regenie/meta.random.out \
  --output-prefix runs/phenotype/regenie/meta

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

> [!IMPORTANT] When additional site methods are implemented, run
> REGENIE, SAIGE, and PLINK results as **separate** meta-analyses —
> don't mix methods in one GWAMA run. All contributing sites must share
> the same phenotype, trait definition, and genome build.

See [`docs/gwama_fixed_random.md`](docs/gwama_fixed_random.md) for installation checks, input requirements, the complete command reference, output formats, result interpretation, testing, and troubleshooting.

### 9. Notes & best practices

- Use **one Brev instance per NVFLARE client**
- Always run the NVFLARE client inside a virtual environment
- Prefer **IAM roles** over static AWS credentials
- Validate GPU availability: `nvidia-smi`
- Use `tmux` or `screen` to keep long-running jobs alive

---

## Synthetic Data Specifications

![Input · Simulation Parameters](https://img.shields.io/badge/INPUT-simulation_parameters-6A4C93?style=for-the-badge)

Generated by `Scripts/generateData.sh` — a full run needs no real cohort data or S3 access. The tables below (marked with the purple badge above) describe what goes *into* the simulation; for what comes *out* of the pipeline, see [Output files](#output-files).

![Samples](https://img.shields.io/badge/samples-88K--110K%2Fsite-6A4C93?style=flat-square)
![SNPs](https://img.shields.io/badge/SNPs-450K--520K%2Fsite-1982C4?style=flat-square)
![Heritability](https://img.shields.io/badge/h²-0.20--0.25-8AC926?style=flat-square)
![Ancestries](https://img.shields.io/badge/ancestries-EUR%20%7C%20EAS%20%7C%20AFR-FFCA3A?style=flat-square)
![Causals](https://img.shields.io/badge/causal%20variants-20%2Fsite-FF595E?style=flat-square)

### Genotypes

| Parameter            | Value                               | Notes                                                        |
|----------------------|-------------------------------------|--------------------------------------------------------------|
| Format               | PLINK binary (`.bed`/`.bim`/`.fam`) | Standard genetic format                                      |
| Variants per site    | 450K–520K SNPs                      | Master grid: 520K; sites subset to their own count           |
| Samples per site     | 88K–110K individuals                | Site 1: 100K, Site 2: 95K, Site 3: 110K                      |
| Chromosomes          | 22 autosomes                        | hg38 assumed (implicit in LDAK)                              |
| Build                | hg38                                | Not explicit in the script; inferred from the PLINK standard |
| MAF range            | 0.01–0.50 (ancestry-dependent)      | EUR: 0.01–0.50; EAS: 0.01–0.45; AFR: 0.005–0.50              |
| LD structure         | Realistic                           | Generated by LDAK; reflects linkage disequilibrium           |
| Population structure | 1–3 subpopulations per site         | Site 1: 1 (homogeneous); Site 2: 2; Site 3: 3                |

### Phenotype

| Parameter         | Value                             | Notes                                                               |
|-------------------|-----------------------------------|---------------------------------------------------------------------|
| File format       | Space-delimited (`FID IID Pheno`) | Standard PLINK phenotype format                                     |
| Trait             | Parkinson's disease (binary)      | 0 = control, 1 = case                                               |
| Prevalence        | 1%                                | Realistic for elderly populations                                   |
| Heritability (h²) | 0.20–0.25 on liability scale      | Site 1: 0.25; Site 2: 0.22; Site 3: 0.20                            |
| Causal variants   | 20 per site                       | 12 shared (identical across sites) + 8 site-specific                |
| Effect size model | LDAK-Thin (power = −0.25)         | Realistic allelic architecture                                      |
| Effect scale      | 1.00, 0.85, 0.70 across sites     | Applied to the first half of shared causals; produces heterogeneity |
| Fixed effects     | Yes (`USE_FIXED_EFFECTS=1`)       | Deterministic, reproducible effects; first half scaled by ancestry  |

### Covariates

| Parameter          | Value                                           | Notes                                                |
|--------------------|-------------------------------------------------|------------------------------------------------------|
| File               | `site{N}_geno.covar`                            | Auto-generated by LDAK `--make-snps`                 |
| Variables          | Age, sex, principal components                  | Exact covariates depend on the LDAK output structure |
| Variance explained | \~10% of phenotypic variation (`COVAR_HER=0.1`) | Controls confounding; realistic level                |

### Shared causal architecture

| Parameter              | Value                                    | Notes                                                                |
|------------------------|------------------------------------------|----------------------------------------------------------------------|
| Shared causal variants | 12                                       | Identical positions across all three sites                           |
| Selection method       | Fixed fractional positions in the genome | Every site generates 520K variants; same physical locations          |
| Effect sizes (shared)  | ±(0.60–0.84)                             | Deterministic function of position; alternating sign                 |
| Ancestry scaling       | 1.00 (EUR), 0.85 (EAS), 0.70 (AFR)       | Applied only to the first 6 shared causals to preserve detectability |

### Site-specific causal architecture

| Parameter                      | Value                                        | Notes                                                  |
|--------------------------------|----------------------------------------------|--------------------------------------------------------|
| Site-specific causals per site | 8                                            | Disjoint between sites; non-overlapping variants       |
| Selection method               | Deterministic band sampling (offset by site) | Spreads variants across the genome; ensures no overlap |
| Effect sizes (site-specific)   | ±(0.70–0.91)                                 | Larger than shared causals; produces local signal      |
| Total causals per site         | 20                                           | 12 shared + 8 specific                                 |

### Key design decisions

1.  **Master grid approach** — all sites generate the same 520K-variant skeleton, then subset to their own count. This puts shared causal variants at identical physical positions without any coordination step. 2.  **Ancestry heterogeneity** — MAF and effect scales differ per site, producing realistic I² (heterogeneity) values in meta-analysis while keeping shared signals detectable. 3.  **Fixed effects** — deterministic and reproducible across runs, which makes validation and debugging tractable. 4.  **Power parameter** — −0.25 reflects realistic allelic architecture, where common variants have smaller effects than rare variants.

---

## Detailed Setup Instructions
### Prerequisites

See [`docs/SYSTEM_REQUIREMENTS.md`](docs/SYSTEM_REQUIREMENTS.md) for the full list of required tools and versions (Python, R, PLINK, REGENIE, GWAMA, LDAK, and OS support), plus
[`docs/requirements.txt`](docs/requirements.txt) for Python package pins.

### REGENIE analysis workflow

`Scripts/local_script_start_gwas.sh` orchestrates the GWAS workflow at a single site: it reads configuration from `~/fedgx_config/config.json`, prepares phenotype and covariate files compatible with REGENIE from raw LDAK output, then runs REGENIE's two-step association analysis. It converts the resulting summary statistics to GWAMA format for aggregation across multiple sites via the NVFLARE federated learning framework. The final output is a standardized summary statistics file that gets streamed back to a central server for meta-analysis.

### Understanding results

#### Output files

![Output · Pipeline Results](https://img.shields.io/badge/OUTPUT-pipeline_results-FF6B35?style=for-the-badge)

These are the files the pipeline *produces*, stage by stage (contrast with the purple-badged [input specification tables](#synthetic-data-specifications) above).

**1. Synthetic data generation — `generateData.sh`**

| File                         | Contents                                   |
|------------------------------|--------------------------------------------|
| `site{N}_geno.bed/.bim/.fam` | PLINK genotypes for the site               |
| `site{N}_geno.covar`         | Covariates (age, sex, PCs)                 |
| `site{N}_pheno.pheno`        | Phenotype (binary trait)                   |
| `site{N}_pheno.liab`         | Underlying liability scores                |
| `site{N}_pheno.effects`      | True simulated effect sizes (ground truth) |
| `site{N}.causals`            | List of causal variant IDs                 |
| `shared_causals.txt`         | Causal variants shared across all sites    |

**2. Per-site REGENIE + standardization — `local_script_start_gwas.sh` → `gwas_cohort_qc_with_gwama.R`**

| File                           | Contents                                                       |
|--------------------------------|----------------------------------------------------------------|
| `{prefix}_pred.list`, `*.loco` | REGENIE step 1 (LOCO model)                                    |
| `{prefix}_{phenotype}.regenie` | REGENIE step 2 raw summary stats                               |
| `{cohort}.cleaned.txt`         | QC'd, harmonized summary stats                                 |
| `{cohort}.GWAMA.txt`           | Standardized GWAMA input (the file that travels to the server) |
| `{cohort}.qc_summary.txt`      | QC pass/fail counts                                            |
| `{cohort}.qc_log.txt`          | Full QC log                                                    |
| `{cohort}.pz_summary.txt`      | P-vs-Z consistency check summary                               |
| `{cohort}.se_n_summary.txt`    | SE-vs-N consistency check summary                              |
| `{cohort}.PZ.png`              | P-Z consistency plot                                           |
| `{cohort}.QQ.png`              | Per-site QQ plot                                               |

**3. Central GWAMA meta-analysis — `run_gwama.sh`**

| File                                      | Contents                                                           |
|-------------------------------------------|--------------------------------------------------------------------|
| `{output_root}.gwama.in`                  | Cohort file list passed to GWAMA                                   |
| `{output_root}.out`                       | Meta-analysis result (single-model runs)                           |
| `{output_root}.fixed.out` / `.random.out` | Fixed/random-effects results (`--model both`)                      |
| `{output_root}.log.out`                   | GWAMA's own run log                                                |
| `{output_root}.comparison.tsv`            | Marker-level FFX vs RFX comparison (via `compare_gwama_models.py`) |

**4. Manhattan plots — `plot_gwama_manhattan.py`**

| File                                   | Contents                                                      |
|----------------------------------------|---------------------------------------------------------------|
| `{output_prefix}.fixed.manhattan.png`  | Fixed-effects Manhattan plot                                  |
| `{output_prefix}.random.manhattan.png` | Random-effects Manhattan plot, same axis scale for comparison |

**5. Interactive dashboard — `build_fedx_dashboard.py` + `fedx_dashboard.html`**

| File           | Contents                                                           |
|----------------|--------------------------------------------------------------------|
| `fedx_data.js` | All site + meta results serialized for the dashboard               |
| `index.html`   | Copy of `fedx_dashboard.html`, only written with `--copy-template` |

`fedgx_meta_aggregator.py` and `serverSide_job.py` don't produce new file types — they orchestrate stages 3–5 in sequence from the server side.

---

## Technologies
| Tool                                                 | Role                                                                                                                                                          |
|------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [LDAK](https://dougspeed.com/) v6.1                  | Data generation                                                                                                                                               |
| [REGENIE](https://rgcgithub.github.io/regenie/) v4.1 | GWAS analysis                                                                                                                                                 |
| [GWAMA](https://genomics.ut.ee/en/tools/gwama)       | Meta-analysis                                                                                                                                                 |
| Docker                                               | Containerization                                                                                                                                              |
| AWS S3                                               | Data storage                                                                                                                                                  |
| NVIDIA FLARE 2.7.1                                   | Federated learning framework                                                                                                                                  |
| `Scripts/client.py` / `Scripts/serverSide_job.py`    | NVFLARE client/server glue: client runs the site's GWAS and standardizes it to GWAMA format; server collects the GWAMA files and hands them to `run_gwama.sh` |
| Brev                                                 | Compute for distributed sites                                                                                                                                 |

---

## References
### Software citations

- **LDAK:** Speed et al. (2020). Improved heritability estimation from genome-wide SNPs. *Nature Genetics*. <https://doi.org/10.1038/s41588-019-0530-8>
- **REGENIE:** Mbatchou et al. (2021). Computationally efficient whole-genome regression for quantitative and binary traits. *Nature Genetics*. <https://doi.org/10.1038/s41588-021-00870-7>
- **GWAMA:** Mägi et al. (2010). GWAMA: software for genome-wide association meta-analysis. *BMC Bioinformatics*. <https://doi.org/10.1186/1471-2105-11-288>
- **This project:** *(add citation when published)*

### Documentation links

- [NVFLARE Documentation](https://nvflare.readthedocs.io/)
- [FedGen Repository](https://github.com/collaborativebioinformatics/FedGen)
- [FedGX Repository](https://github.com/collaborativebioinformatics/FedGX)
- [Brev Platform](https://brev.dev)
- [REGENIE Documentation](https://rgcgithub.github.io/regenie/)
- [LDAK Documentation](https://dougspeed.com/)
- [PLINK File Formats](https://www.cog-genomics.org/plink/1.9/formats)

---

## License

Data and scripts: MIT License (see repository root).

---

## Team
| Contributor        | Sections | Affiliation                                                                                                                                                                                                              |
|--------------------|----------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Marlene Rietz      | 1–3      | Steno Diabetes Center Odense, Odense, Denmark; P1 Pioneer Center for Artificial Intelligence, University of Copenhagen, Copenhagen, Denmark; Department of Laboratory Medicine, Karolinska Institutet, Stockholm, Sweden |
| Allan Lind-Thomsen | 4        | OPEN, Odense University Hospital, Denmark                                                                                                                                                                                |
| Xiaoping Wu        | 5        | Department of Obstetrics and Gynecology, Institute of Clinical Sciences, Sahlgrenska Academy, University of Gothenburg, Gothenburg, Sweden                                                                               |
| Moh Sallam         | 6        | Center for Quantitative Genetics and Genomics and Pioneer Center for SMARTbiomed, Aarhus University, Denmark                                                                                                                      |
| Pravesh Parekh     | 7–8      | J. Craig Venter Institute, San Diego, California, USA; Centre for Precision Psychiatry, University of Oslo, Oslo, Norway                                                                                                 |
with thanks to Ziyue Xu and Holger Roth from NVIDIA!
---

<p align="center">

<sub>Built at the Nordic Biobank × NVIDIA Hackathon</sub>

</p>
