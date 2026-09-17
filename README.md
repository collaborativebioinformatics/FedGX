# FEDGX: Federated Learning Software for Multi-Tool Genome-Wide Association Studies across Cohort Sites
## Product of the Nordic Biobank x NVIDIA Hackathon

  ![Logo](docs/logo.jpeg)

  ---

## Table of Contents
- [Summary](#summary)
- [Flowchart](#flowchart)
- [Quickstart](#quickstart)
- [Data Documentation](#data-specifications)
- [Detailed Setup](#detailed-setup-instructions)
- [References](#references)

---

# Summary

### At the Nordic Biobank x NVIDIA Hackathon, we aimed to develop software for federated GWAS across three sites. 
For development, we tested this approach in the HUNT Cloud and across two BREV sites. 

We extended the currently available code from the FedGen repository. 
https://github.com/collaborativebioinformatics/FedGen

____________________________

# Flowchart
  ![Federated GWAS architecture](docs/FederationFigure_MR.png)
____________________________

# Quickstart

## 1. Start NVFLARE Dashboard and FL Server

## 2. Start NVFLARE Client on Brev and HUNT Cloud

### 2.1 Create GPU Instance on Brev

On the **Brev website**:

* Create **1 GPU instance** per site
* Example configuration:
  * Name: `site1`
  * GPU: **1× NVIDIA L4**
  * CPU: **16 cores**
  * RAM: **64 GB**

### 2.2 Connect to the Instance

```bash
brev shell site1
```

Use terminal multiplexer to ensure connection persistence (Optional but recommended)

```bash
tmux new -s nvflare
```

### 2.3 Python Environment Setup

```bash
python3 -m venv venv_nvflare
source venv_nvflare/bin/activate

pip install nvflare[PT] torch torchvision tensorboard
```

Verify installation:

```bash
nvflare --version
```

## 3. Copy and Start NVFLARE Client Startup Kit

### 3.1 Copy Client Kit from Local Machine

On **local machine**:

```bash
brev copy <local_path_to_client_kit> site1:<remote_path>
```

On **Brev instance**:

```bash
sudo apt update
sudo apt install -y unzip

unzip -d <client_name> -P <PIN> <client_kit.zip>
cd <client_name>
```

### 3.2 Start NVFLARE Client

```bash
./startup/start.sh
```

Check logs to confirm successful connection to the NVFLARE server/dashboard.


## 4. Clone FedGX Repository

```bash
git clone https://github.com/collaborativebioinformatics/FedGX
chmod +x FedGen/scripts/*.sh
```


## 5. Download Site Data from S3

[SPECIFY METHOD HERE; SYNTHETIC DATA CODE AVAILABLE]


## 6. Run Regenie Per Site (Outside NVFLARE)

Run Regenie independently per site (not through NVFLARE) to verify all dependencies are working:

```bash
cd ~/data
./../FedGen/scripts/run_regenie_site.sh <siteNumber>
```

Monitor logs and outputs to confirm successful completion.

**Runtime:** ~30-45 minutes total
- Step 1 (LOCO model): 15-30 min
- Step 2 (association testing): 10-20 min


## 7. Run Federated GWAS Job (NVFLARE)

Instead of running REGENIE independently on each site and manually aggregating results, you can submit a federated GWAS job that automates the entire workflow across all sites using NVIDIA FLARE.

The federated job handles:
- Distributing analysis scripts to all clients
- Running local GWAS analysis using REGENIE on each site
- Collecting summary statistics from all sites
- Performing meta-analysis using GWAMA on the server

**For complete instructions on submitting federated GWAS jobs, see [`jobs/fed_gwas/README.md`](jobs/fed_gwas/README.md).**


## 8. Run central GWAMA meta-analysis

The central runner supports fixed-effect (`fixed`), random-effects (`random`),
or paired (`both`) GWAMA analyses for binary and quantitative traits. Paired
runs also produce a marker-level FFX/RFX comparison. A separate script creates
matched fixed- and random-effects Manhattan plots.

Before the central run, each site uses
`Scripts/gwas_cohort_qc_with_gwama.R` to QC its REGENIE, SAIGE, or PLINK
summary results and produce `<site>.GWAMA.txt.gz`. The command requires an
explicit mapping from the source columns to chromosome, position, alleles,
frequency, effect, standard error, sample size, and P-value. See the detailed
guide for a complete REGENIE example and requirements for other formats.

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

Run REGENIE, SAIGE, and PLINK results as separate meta-analyses; do not mix
methods in one GWAMA run. All contributing sites must use the same phenotype,
trait definition, and genome build.

See [Central GWAMA fixed/random analysis](docs/gwama_fixed_random.md) for
installation checks, input requirements, complete commands, outputs,
interpretation, testing, and troubleshooting.

## 9. Notes & Best Practices

* Use **one Brev instance per NVFLARE client**
* Always run NVFLARE client inside a virtual environment
* Prefer **IAM roles** over static AWS credentials
* Validate GPU availability:

  ```bash
  nvidia-smi
  ```
* Use `tmux` or `screen` to keep long‑running jobs alive

---
# Synthetic Data Specifications

## Genotypes

| Parameter | Value | Notes |
|-----------|-------|-------|
| Format | PLINK binary (.bed/.bim/.fam) | Standard genetic format |
| Variants per site | 450K–520K SNPs | Master grid: 520K; sites subset to their own count |
| Samples per site | 88K–110K individuals | Site 1: 100K, Site 2: 95K, Site 3: 110K |
| Chromosomes | 22 autosomes | hg38 assumed (implicit in LDAK) |
| Build | hg38 | Not explicit in script; inferred from PLINK standard |
| MAF range | 0.01–0.50 (ancestry-dependent) | EUR: 0.01–0.50; EAS: 0.01–0.45; AFR: 0.005–0.50 |
| LD structure | Realistic | Generated by LDAK; reflects linkage disequilibrium |
| Population structure | 1–3 subpopulations per site | Site 1: 1 (homogeneous); Site 2: 2; Site 3: 3 |

## Phenotype

| Parameter | Value | Notes |
|-----------|-------|-------|
| File format | Space-delimited (FID IID Pheno) | Standard PLINK phenotype format |
| Trait | Parkinson's disease (binary) | 0 = control, 1 = case |
| Prevalence | 1% | Realistic for elderly populations |
| Heritability (h²) | 0.20–0.25 on liability scale | Site 1: 0.25; Site 2: 0.22; Site 3: 0.20 |
| Causal variants | 20 per site | 12 shared (identical across sites) + 8 site-specific |
| Effect size model | LDAK-Thin (power = –0.25) | Realistic allelic architecture |
| Effect scale | 1.00, 0.85, 0.70 across sites | Applied to first half of shared causals; produces heterogeneity |
| Fixed effects | Yes (USE_FIXED_EFFECTS=1) | Deterministic, reproducible effects; first half scaled by ancestry |

## Covariates

| Parameter | Value | Notes |
|-----------|-------|-------|
| File | `site{N}_geno.covar` | Auto-generated by LDAK --make-snps |
| Variables | Age, sex, principal components | Specific covariates depend on LDAK output structure |
| Variance explained | ~10% of phenotypic variation (COVAR_HER=0.1) | Controls confounding; realistic level |

## Shared Causal Architecture

| Parameter | Value | Notes |
|-----------|-------|-------|
| Shared causal variants | 12 | Identical positions across all three sites |
| Selection method | Fixed fractional positions in genome | Every site generates 520K variants; same physical locations |
| Effect sizes (shared) | ±(0.60–0.84) | Deterministic function of position; alternating sign |
| Ancestry scaling | 1.00 (EUR), 0.85 (EAS), 0.70 (AFR) | Applied only to first 6 shared causals to preserve detectability |

## Site-Specific Causal Architecture

| Parameter | Value | Notes |
|-----------|-------|-------|
| Site-specific causals per site | 8 | Disjoint between sites; non-overlapping variants |
| Selection method | Deterministic band sampling (offset by site) | Spreads variants across genome; ensures non-overlap |
| Effect sizes (site-specific) | ±(0.70–0.91) | Larger than shared causals; produces local signal |
| Total causals per site | 20 | 12 shared + 8 specific |

## Key Design Decisions

1. **Master Grid Approach**: All sites generate the same 520K variant skeleton, then subset to their own count. Ensures shared causal variants occupy identical physical positions without coordination.

2. **Ancestry Heterogeneity**: MAF and effect scales differ per site, producing realistic I² values (heterogeneity) in meta-analysis while keeping shared signals detectable.

3. **Fixed Effects**: Deterministic, reproducible across runs; enables validation and debugging.

4. **Power Parameter**: –0.25 reflects realistic allelic architecture where common variants have smaller effects than rare variants.

---


# Detailed Setup Instructions

## Prerequisites


---

## Download Workflow

```bash
# 1. Clone repository (if not already done)
git clone https://github.com/collaborativebioinformatics/FedGen.git
cd FedGen

# 2. Download your assigned site (e.g., Site 3)
./scripts/download_site_from_s3.sh 3

# 3. Verify download
ls -lh data/simulated_sites/site3/
# Should show ~15 GB total:
# - site3_geno.bed (~12-13 GB)
# - site3_geno.bim (~10-20 MB)
# - site3_geno.fam (~2-3 MB)
# - site3_pheno.pheno (~2 MB)
# - site3_geno.covar (~5-10 MB)
```

---


## REGENIE Analysis Workflow

ALLAN?



---

### Run Analysis


---

## Understanding Results#

XIAOPING?

### Output Files

### Association Results Format
INSERT

**Key columns:**
INSERT


### Find Genome-Wide Significant Hits

### Manhattan Plot (R)



---




# Technologies

- **Data Generation:** [LDAK](https://dougspeed.com/) v6.1
- **GWAS Analysis:** [REGENIE](https://rgcgithub.github.io/regenie/) v4.1
- **Meta-Analysis:** [GWAMA](https://genomics.ut.ee/en/tools/gwama)
- **Containerization:** Docker
- **Data Storage:** AWS S3
- **FL Framework:** NVIDIA FLARE 2.7.1
- **Compute:** Brev instances for distributed sites

---

# References

## Software Citations

- **LDAK:** Speed et al. (2020). Improved heritability estimation from genome-wide SNPs. *Nature Genetics*. https://doi.org/10.1038/s41588-019-0530-8

- **REGENIE:** Mbatchou et al. (2021). Computationally efficient whole-genome regression for quantitative and binary traits. *Nature Genetics*. https://doi.org/10.1038/s41588-021-00870-7

- **GWAMA:** Mägi et al. (2010). GWAMA: software for genome-wide association meta-analysis. *BMC Bioinformatics*. https://doi.org/10.1186/1471-2105-11-288

- **This project:** [Add citation when published]

## Documentation Links

- **NVFLARE Documentation:** [https://nvflare.readthedocs.io/](https://nvflare.readthedocs.io/)
- **FedGen Repository:** [https://github.com/collaborativebioinformatics/FedGen](https://github.com/collaborativebioinformatics/FedGen)
- **Brev Platform:** [https://brev.dev](https://brev.dev)
- **REGENIE Documentation:** https://rgcgithub.github.io/regenie/
- **LDAK Documentation:** https://dougspeed.com/
- **PLINK File Formats:** https://www.cog-genomics.org/plink/1.9/formats

---

---

# License

Data and scripts: MIT License (see repository root)

# Contributors

---

# Team

Marlene Rietz (1-3)
Allan Lind-Thomsen (4)
Xiaoping Wu (5)
Moh Sallam
(6)
Pravesh Parekh (7-8)

## Affiliations

1.  Steno Diabetes Center Odense, Odense, Denmark

2.  P1 Pioneer Center for Artificial Intelligence, University of
    Copenhagen, Copenhagen Denmark

3.  Department of Laboratory Medicine, Karolinska Institutet, Stockholm,
    Sweden

4.  [INSERT ALLAN]

5.  Department of Obstetrics and Gynecology,Institute of Clinical
    Sciences, Sahlgrenska Academy, University of Gothenburg, Gothenburg,
    Sweden

6.  Center for Quantitative Genetics and Genomics and Pionner Center for
    Smartbiomed, Aarhus University

7.  J. Craig Venter Institute, San Diego, California, USA

8.  Centre for Precision Psychiatry, University of Oslo, Oslo, Norway




