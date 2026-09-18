# FedGX System & Tool Requirements

This document specifies the system-level tools, compiled binaries, and their required versions for the FedGX federated GWAS pipeline.

## Operating System

- **Linux** (Ubuntu 20.04+, Debian 10+, or compatible)
- **macOS** (10.15+) — note: GPU support is limited
- **HPC/Cloud**: HUNT Cloud, Brev GPU instances

## Programming Languages & Runtimes

| Tool     | Version | Purpose                           | Install                                     |
|----------|---------|-----------------------------------|---------------------------------------------|
| Python   | 3.10+   | Federated learning, analysis      | `apt install python3.10`                    |
| R        | 4.0+    | GWAMA output processing, plotting | `apt install r-base`                        |
| Bash     | 4.0+    | Shell scripts for workflows       | Standard on most Linux systems               |

## Build & Development Tools

| Tool      | Version | Purpose                                | Install                              |
|-----------|---------|----------------------------------------|--------------------------------------|
| gcc/g++   | 9.0+    | Compile Python packages & binaries     | `apt install build-essential`        |
| Make      | 4.0+    | Build system for some dependencies     | Included in `build-essential`         |
| git       | 2.25+   | Clone FedGX repository                 | `apt install git`                    |

## Genomics Tools (Compiled Binaries)

Install these **outside** the Python virtual environment or add to PATH.

### REGENIE v4.1
- **Purpose**: Whole-genome regression for GWAS at individual sites
- **Download**: https://github.com/rgcgithub/regenie/releases/tag/v4.1.0.1
- **Install**:
  ```bash
  wget https://github.com/rgcgithub/regenie/releases/download/v4.1.0.1/regenie_v4.1.0.1.gz_x86_64_Linux.tar.gz
  tar -xzf regenie_v4.1.0.1.gz_x86_64_Linux.tar.gz
  sudo mv regenie /usr/local/bin/
  regenie --help  # verify
  ```
- **Runtime**: 15–45 minutes per site (depends on sample size and SNP count)
- **Dependencies**: None (precompiled)

### GWAMA (Genome-Wide Association Meta-Analysis)
- **Purpose**: Fixed/random-effects meta-analysis across sites
- **Download**: https://genomics.ut.ee/en/tools/gwama
- **Install**:
  ```bash
  wget https://genomics.ut.ee/en/tools/gwama
  chmod +x gwama
  sudo mv gwama /usr/local/bin/
  gwama --help  # verify
  ```
- **Runtime**: 1–10 minutes for 3-site meta-analysis
- **Language**: Fortran (precompiled binary provided)

### LDAK v6.1
- **Purpose**: Simulating genotypes and LD structure for test data
- **Download**: https://dougspeed.com/
- **Install**:
  ```bash
  # Follow download instructions from LDAK website (requires registration)
  tar -xzf ldak5.linux.tar.gz
  sudo mv ldak5.linux /usr/local/bin/ldak
  ldak --help  # verify
  ```
- **Runtime**: 30–60 minutes for generating 3 sites with ~500K SNPs each
- **Use**: `Scripts/generateData.sh` calls LDAK for data generation

### PLINK / PLINK2
- **Purpose**: Genetic data format conversion and QC
- **Download**: 
  - PLINK v1.9: https://www.cog-genomics.org/plink/1.9/
  - PLINK2: https://www.cog-genomics.org/plink/2.0/
- **Install**:
  ```bash
  # PLINK v1.9
  wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230616.zip
  unzip plink_linux_x86_64_20230616.zip
  sudo mv plink /usr/local/bin/
  plink --version  # verify
  ```
- **Purpose in FedGX**: Data validation, format checks
- **Optional**: Required only if converting to/from non-PLINK formats

## Session Management

| Tool      | Version | Purpose                                    | Install/Setup                       |
|-----------|---------|--------------------------------------------|------------------------------------|
| tmux      | 3.0+    | Session multiplexer for persistent jobs    | `apt install tmux`                  |
| screen    | 4.0+    | Alternative session manager                | `apt install screen`                |

## GPU & CUDA (for Brev instances & GPU-enabled HPC)

| Tool              | Version | Purpose                           | Install                            |
|-------------------|---------|-----------------------------------|------------------------------------|
| CUDA Toolkit      | 11.8+   | GPU compute runtime               | https://developer.nvidia.com/cuda-11-8-0-download-archive |
| cuDNN             | 8.0+    | Deep learning GPU library         | https://developer.nvidia.com/cudnn |
| nvidia-utils      | Latest  | GPU driver utilities              | `apt install nvidia-utils`         |
| nvidia-driver     | 525+    | NVIDIA GPU driver                 | `apt install nvidia-driver-525`    |

**Verify GPU setup**:
```bash
nvidia-smi  # should show GPU memory and driver version
```

## Python Virtual Environment Setup

Create an isolated environment for FedGX:

```bash
# 1. Create virtual environment
python3 -m venv venv_nvflare

# 2. Activate it
source venv_nvflare/bin/activate

# 3. Upgrade pip, setuptools, wheel
pip install --upgrade pip setuptools wheel

# 4. Install Python dependencies from requirements.txt
pip install -r requirements.txt

# 5. Verify key installations
nvflare --version
python3 -c "import torch; print(torch.__version__)"
python3 -c "import regenie" 2>&1 || echo "REGENIE must be in PATH or installed separately"
```

## Typical Ubuntu Installation Script

```bash
#!/bin/bash
set -e

echo "Installing FedGX dependencies on Ubuntu..."

# System updates
sudo apt update
sudo apt upgrade -y

# Build tools & dev libraries
sudo apt install -y build-essential python3-dev python3-venv git

# Python 3.10+ (if not already installed)
sudo apt install -y python3.10 python3.10-venv python3.10-dev

# R (for GWAMA output processing)
sudo apt install -y r-base r-base-dev

# Session managers
sudo apt install -y tmux screen

# GPU support (if applicable)
# sudo apt install -y nvidia-driver-525 nvidia-utils

echo "System dependencies installed. Now install Python packages:"
echo "  python3 -m venv venv_nvflare"
echo "  source venv_nvflare/bin/activate"
echo "  pip install -r requirements.txt"
echo ""
echo "Then manually download and install:"
echo "  - REGENIE v4.1"
echo "  - GWAMA"
echo "  - LDAK v6.1 (optional, for synthetic data generation)"
echo "  - PLINK (optional, for format conversion)"
```

## Verification Checklist

Run these commands to verify all dependencies are installed:

```bash
# Python & virtual environment
python3 --version
python3 -m venv --help

# Build tools
gcc --version
g++ --version
git --version

# AWS CLI
aws --version

# After activating venv_nvflare:
source venv_nvflare/bin/activate
python3 -c "import nvflare; print('NVFLARE OK')"
python3 -c "import torch; print('PyTorch OK')"
python3 -c "import pandas; print('Pandas OK')"
python3 -c "import numpy; print('NumPy OK')"

# Genomics tools (should be in PATH)
regenie --help
gwama --help
ldak --help  # if installed
plink --version  # if installed

# Session managers
tmux --version
screen --version

# GPU (if applicable)
nvidia-smi
```

## Common Installation Issues & Solutions

### Issue: `pip install nvflare[PT]` fails on GPU instance
**Solution**: Ensure CUDA 11.8+ and cuDNN 8.0+ are installed before installing PyTorch.
```bash
nvcc --version  # verify CUDA
```

### Issue: `regenie` not found when running scripts
**Solution**: Add REGENIE to PATH or provide full path in scripts.
```bash
export PATH=$PATH:/usr/local/bin  # if regenie is there
# Or edit Scripts/run_regenie_site.sh to use full path
```

### Issue: Out of memory on Brev instance
**Solution**: Ensure instance size matches data. For 100K samples × 500K SNPs, use at least:
- **16 GB RAM** (minimum)
- **64 GB RAM** (recommended for parallel processing)
- **1× NVIDIA L4 GPU** (optional, for acceleration)



## References

- REGENIE: https://rgcgithub.github.io/regenie/
- GWAMA: https://genomics.ut.ee/en/tools/gwama
- LDAK: https://dougspeed.com/
- PLINK: https://www.cog-genomics.org/plink/
- NVIDIA FLARE: https://nvflare.readthedocs.io/
- Brev platform: https://brev.dev
- HUNT Cloud: https://www.ntnu.edu/mh/hunt/hunt-cloud

---

*Last updated: September 2026*
*For FedGX at Nordic Biobank × NVIDIA Hackathon*
