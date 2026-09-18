---
marp: true
theme: default
paginate: true
---

# FedGX
### Federated GWAS Across Cohort Sites

Nordic Biobank × NVIDIA Hackathon

![width:180px](docs/logo.jpeg)

---

## The problem

- Genome-wide association studies need **large, diverse cohorts** to find real signal
- But cohort data usually can't leave the site it came from — privacy, governance, legal restrictions
- Result: most GWAS run on a single biobank, missing ancestries it doesn't have

**Question:** can we analyze across sites without centralizing the raw genotypes?

---

## What we built

**FedGX** — a federated GWAS pipeline where each site:

1. Keeps its genotypes local, always
2. Runs its own association analysis (REGENIE)
3. Sends only **summary statistics** to a central server
4. Gets combined via meta-analysis (GWAMA)

Built on [NVIDIA FLARE](https://nvidia.github.io/NVFlare/), extending the [FedGen](https://github.com/collaborativebioinformatics/FedGen) hackathon scaffolding.

---

## Architecture

![width:900px](docs/FederationFigure_MR.png)

3 sites, 1 ancestry proxy each (EUR / EAS / AFR) — nothing but summary statistics ever leaves a site.

---

## Pipeline, stage by stage

| Stage | What happens |
|---|---|
| 1. Synthetic data | LDAK generates 3 sites' genotypes + phenotypes locally |
| 2. Per-site GWAS | REGENIE step 1 + 2, then QC/standardize to GWAMA format |
| 3. Federation | NVFLARE client sends the standardized file to the server |
| 4. Meta-analysis | GWAMA combines all sites — fixed **and** random effects |
| 5. Visualization | Matched Manhattan plots + an interactive dashboard |

---

## Synthetic data: not just random noise

To make the meta-analysis meaningful, the 3 sites share a **master SNP grid** and:

- **12 causal variants shared** across all 3 sites (with ancestry-scaled effect sizes)
- **8 causal variants unique** to each site
- Different MAF spectra per site (EUR/EAS/AFR proxy)
- Different heritability, sample size, and population substructure per site

This produces realistic **heterogeneity** — the same signal, weakened or shifted by ancestry — instead of a toy case where every site is identical.

---

## Fixed vs. random effects — and why both

- **Fixed-effect (FFX):** assumes every site estimates the *same* true effect — most power, but wrong if sites genuinely differ
- **Random-effects (RFX):** allows the true effect to vary by site — more honest under heterogeneity, less power

FedGX runs **both** by default and produces a marker-level comparison automatically:

- Where do FFX and RFX estimates diverge?
- Does that line up with high I² (heterogeneity)?

---

## Results & visualization

- Matched **Manhattan plots** for fixed- and random-effects results, same axis scale, side by side
- An **interactive dashboard** (`fedx_dashboard.html`) built from real pipeline output — no server needed, just open it in a browser
- Full per-site QC trail: P-Z consistency, SE-vs-N consistency, QQ plots

<video src="docs/dashboard_demo.mov" controls width="700"></video>

*Live walkthrough of the interactive dashboard (`fedx_dashboard.html`) built from real pipeline output.*

---

## Technology stack

| | |
|---|---|
| **Data generation** | LDAK |
| **Per-site GWAS** | REGENIE |
| **Meta-analysis** | GWAMA |
| **Federated learning** | NVIDIA FLARE |
| **Compute** | Brev, HUNT Cloud |
| **Storage** | AWS S3 |

---

## What we delivered vs. the original plan

| Planned | Status |
|---|---|
| Synthetic 3-site multi-ancestry data | ✅ Done |
| Per-site REGENIE client pipeline | ✅ Done |
| Central summary-stat standardization | ✅ Done |
| GWAMA meta-analysis — fixed effects | ✅ Done |
| GWAMA meta-analysis — random effects | ✅ Done (was flagged buggy at planning time) |
| Visualization (originally "bonus") | ✅ Done — plots + interactive dashboard |
| PLINK / GCTA / SAIGE site drivers | 🔜 Scaffolded, not yet implemented |
| Polygenic risk scores (PRS) | 🔜 Open question, out of scope this round |

---

## Try it yourself

```bash
git clone https://github.com/collaborativebioinformatics/FedGX
cd FedGX
./Scripts/generateData.sh 1   # no real cohort data needed
```

Full setup: [README.md](README.md)
Fixed/random-effects deep dive: [docs/gwama_fixed_random.md](docs/gwama_fixed_random.md)

---

## Team

| | |
|---|---|
| **Marlene Rietz** | Steno Diabetes Center Odense · UCPH · Karolinska Institutet |
| **Allan Lind-Thomsen** | OPEN, Odense University Hospital |
| **Xiaoping Wu** | University of Gothenburg |
| **Moh Sallam** | Aarhus University |
| **Pravesh Parekh** | J. Craig Venter Institute · University of Oslo |

Built at the Nordic Biobank × NVIDIA Hackathon

---

# Thank you

**github.com/collaborativebioinformatics/FedGX**
