# FEDGX: Federated Learning Software for Multi-Tool Genome-Wide Association Studies across Cohort Sites: Product of the Nordic Biobank x NVIDIA Hackathon

# Contributors

Marlene Rietz (1-3), Allan Lind-Thomsen (4), Xiaoping Wu (5), Moh Sallam
(6), Pravesh Parekh (7-8)

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


# Summary

At the Nordic Biobank x NVIDIA Hackathon, we aim to develop software for federated GWAS across three sites. 
For development, we test this approach in the HUNT Cloud and across two BREV sites. 

We will extend currently available code from the FedGen repository. 
https://github.com/collaborativebioinformatics/FedGen

____________________________

# Flowchart
  ![Federated GWAS architecture](docs/FederationFigure_MR.png)
____________________________
 

____________________________

## Initial Plan 

1. Synthetic dataset creation across 3 sites (genotype and phenotype)

- Use LDAK from DougSpeed.com
- Generate genotype across three different populations (ancestries)
- Genotype/phenotype mapping

2. Client pipeline for GWAS software in *Site 1, 2, 3*

- Extend to handle PLINK, GCTA, SAIGE, and custom approaches, over and above REGENIE
- Edits to server-side GWAS code to handle different GWAS calls
- Possibly LD structure handling

3. Standardization of summary stats in *Central Analytical Engine*

4. Meta-analyses in *Central Analytical Engine* using GWAMA

- FFX — currently implemented
- RFX — currently buggy
- Optionally develop LD structure weighted meta-analyses

5. Visualisation component (bonus, not core scope)

6. PRS (open question — not yet scoped)

____________________________

---

## Table of Contents
- [Quickstart](#quickstart----server-and-clients-configuration)
- [Data Documentation](#data-specifications)
- [Detailed Setup](#detailed-setup-instructions)
- [Troubleshooting](#troubleshooting)
- [References](#references)

---

____________________________

# Quickstart -- Server and Clients Configuration

## 1. Start NVFLARE Dashboard and FL Server on AWS



---

## 2. Start NVFLARE Client on Brev and HUNT Cloud

### 2.1 Create GPU Instance on Brev

On the **Brev website**:

* Create **1 GPU instance** per site
* Example configuration:
  * Name: `site1`
  * GPU: **1× NVIDIA L4**
  * CPU: **16 cores**
  * RAM: **64 GB**

---

### 2.2 Connect to the Instance

```bash
brev shell site1
```

Use terminal multiplexer to ensure connection persistence (Optional but recommended)

```bash
tmux new -s nvflare
```

---

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

---

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

---

## 4. Install AWS CLI on Each Brev Instance

From your **home directory**:

```bash
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
```

Verify:

```bash
aws --version
```

---

### 4.1 Configure AWS Credentials (Securely)

```bash
aws configure
```

Use **one** of the following secure approaches:

* IAM role attached to the instance (**recommended**)
* Environment variables (`AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`)
* AWS credentials file

Example (DO NOT hardcode secrets):

```
AWS Access Key ID:     <YOUR_ACCESS_KEY>
AWS Secret Access Key: <YOUR_SECRET_KEY>
Default region name:  None
Default output format: None
```

---

## 5. Clone FedGX Repository

```bash
git clone https://github.com/collaborativebioinformatics/FedGX
chmod +x FedGen/scripts/*.sh
```

---

## 6. Download Site Data from S3

[SPECIFY METHOD HERE; SYNTHETIC DATA CODE AVAILABLE]











