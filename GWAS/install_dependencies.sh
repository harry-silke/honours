#!/bin/bash
#
# Dependency Installation Script for Kelp GWAS Pipeline
#
# This script installs all required software dependencies for the GWAS pipeline
# Supports: Ubuntu/Debian, CentOS/RHEL, and macOS
#
# Usage:
#   bash install_dependencies.sh [--method conda|manual]
#
# Author: GWAS Pipeline
# Date: 2026-01-26

set -e  # Exit on error

echo "============================================================================="
echo "Kelp GWAS Pipeline - Dependency Installation"
echo "============================================================================="
echo

# Parse arguments
INSTALL_METHOD="conda"  # Default to conda

while [[ $# -gt 0 ]]; do
  case $1 in
    --method)
      INSTALL_METHOD="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: $0 [--method conda|manual]"
      exit 1
      ;;
  esac
done

echo "Installation method: $INSTALL_METHOD"
echo

# =============================================================================
# Function Definitions
# =============================================================================

check_command() {
  if command -v "$1" &> /dev/null; then
    echo "✓ $1 is installed"
    return 0
  else
    echo "✗ $1 is not installed"
    return 1
  fi
}

# =============================================================================
# Check Existing Installations
# =============================================================================

echo "--- Checking existing installations ---"
echo

TOOLS=("python3" "R" "fastqc" "fastp" "hisat2" "samtools" "bcftools" "gatk" "plink" "bedtools")

for tool in "${TOOLS[@]}"; do
  check_command "$tool" || true
done

echo

# =============================================================================
# Conda Installation (Recommended)
# =============================================================================

if [ "$INSTALL_METHOD" == "conda" ]; then
  echo "============================================================================="
  echo "Installing via Conda (Recommended)"
  echo "============================================================================="
  echo

  # Check if conda is installed
  if ! command -v conda &> /dev/null; then
    echo "Conda not found. Please install Miniconda or Anaconda first:"
    echo "https://docs.conda.io/en/latest/miniconda.html"
    exit 1
  fi

  # Create conda environment
  ENV_NAME="kelp_gwas"

  echo "Creating conda environment: $ENV_NAME"

  # Check if environment exists
  if conda env list | grep -q "^$ENV_NAME "; then
    echo "Environment $ENV_NAME already exists."
    read -p "Do you want to remove it and create a fresh one? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
      conda env remove -n $ENV_NAME -y
    else
      echo "Updating existing environment..."
    fi
  fi

  # Create/update environment
  echo "Installing packages (this may take several minutes)..."

  conda create -n $ENV_NAME python=3.9 -y

  # Activate environment
  eval "$(conda shell.bash hook)"
  conda activate $ENV_NAME

  # Install bioinformatics tools
  conda install -c bioconda -c conda-forge -y \\
    fastqc \\
    fastp \\
    hisat2 \\
    samtools \\
    bcftools \\
    gatk4 \\
    plink \\
    bedtools \\
    cyvcf2 \\
    pysam \\
    biopython \\
    pybedtools

  # Install Python packages
  conda install -c conda-forge -y \\
    numpy \\
    pandas \\
    scipy \\
    matplotlib \\
    seaborn \\
    pyyaml \\
    tqdm

  pip install statsmodels

  # Install R and R packages
  conda install -c conda-forge -y \\
    r-base \\
    r-essentials \\
    r-data.table \\
    r-dplyr \\
    r-ggplot2 \\
    r-glmnet \\
    r-proc \\
    r-caret \\
    r-survival

  # Install Bioconductor packages in R
  echo "Installing R packages from Bioconductor..."
  Rscript install_R_packages.R

  echo
  echo "============================================================================="
  echo "Installation Complete!"
  echo "============================================================================="
  echo
  echo "To use the pipeline, activate the environment:"
  echo "  conda activate $ENV_NAME"
  echo
  echo "To deactivate:"
  echo "  conda deactivate"
  echo
  echo "To verify installation:"
  echo "  bash install_dependencies.sh --verify"
  echo

# =============================================================================
# Manual Installation
# =============================================================================

elif [ "$INSTALL_METHOD" == "manual" ]; then
  echo "============================================================================="
  echo "Manual Installation Guide"
  echo "============================================================================="
  echo
  echo "Please install the following tools manually:"
  echo
  echo "1. Python 3.8+ with packages:"
  echo "   pip install -r requirements.txt"
  echo
  echo "2. R 4.0+ with packages:"
  echo "   Rscript install_R_packages.R"
  echo
  echo "3. Bioinformatics tools:"
  echo
  echo "   FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"
  echo "   fastp: https://github.com/OpenGene/fastp"
  echo "   HISAT2: http://daehwankimlab.github.io/hisat2/"
  echo "   SAMtools: http://www.htslib.org/"
  echo "   BCFtools: http://www.htslib.org/"
  echo "   GATK: https://gatk.broadinstitute.org/"
  echo "   PLINK: https://www.cog-genomics.org/plink/2.0/"
  echo "   BEDtools: https://bedtools.readthedocs.io/"
  echo
  echo "On Ubuntu/Debian:"
  echo "  sudo apt-get update"
  echo "  sudo apt-get install fastqc samtools bcftools bedtools"
  echo
  echo "On CentOS/RHEL:"
  echo "  sudo yum install samtools bcftools bedtools"
  echo
  echo "On macOS (with Homebrew):"
  echo "  brew install fastqc samtools bcftools bedtools"
  echo
  echo "For GATK and PLINK, download from their respective websites."
  echo

# =============================================================================
# Verification
# =============================================================================

elif [ "$INSTALL_METHOD" == "verify" ]; then
  echo "============================================================================="
  echo "Verifying Installation"
  echo "============================================================================="
  echo

  ALL_INSTALLED=true

  echo "--- Checking Python packages ---"
  python3 -c "import numpy; print('✓ numpy')" || { echo "✗ numpy"; ALL_INSTALLED=false; }
  python3 -c "import pandas; print('✓ pandas')" || { echo "✗ pandas"; ALL_INSTALLED=false; }
  python3 -c "import cyvcf2; print('✓ cyvcf2')" || { echo "✗ cyvcf2"; ALL_INSTALLED=false; }
  python3 -c "import matplotlib; print('✓ matplotlib')" || { echo "✗ matplotlib"; ALL_INSTALLED=false; }

  echo
  echo "--- Checking R packages ---"
  Rscript -e "library(data.table); cat('✓ data.table\n')" 2>/dev/null || { echo "✗ data.table"; ALL_INSTALLED=false; }
  Rscript -e "library(ggplot2); cat('✓ ggplot2\n')" 2>/dev/null || { echo "✗ ggplot2"; ALL_INSTALLED=false; }
  Rscript -e "library(SNPRelate); cat('✓ SNPRelate\n')" 2>/dev/null || { echo "✗ SNPRelate"; ALL_INSTALLED=false; }
  Rscript -e "library(qqman); cat('✓ qqman\n')" 2>/dev/null || { echo "✗ qqman"; ALL_INSTALLED=false; }

  echo
  echo "--- Checking bioinformatics tools ---"
  for tool in "${TOOLS[@]}"; do
    check_command "$tool" || ALL_INSTALLED=false
  done

  echo
  if [ "$ALL_INSTALLED" = true ]; then
    echo "============================================================================="
    echo "All dependencies are installed! ✓"
    echo "============================================================================="
    exit 0
  else
    echo "============================================================================="
    echo "Some dependencies are missing. Please install them."
    echo "============================================================================="
    exit 1
  fi

else
  echo "Unknown installation method: $INSTALL_METHOD"
  echo "Use: conda, manual, or verify"
  exit 1
fi

# =============================================================================
# Create activation script
# =============================================================================

if [ "$INSTALL_METHOD" == "conda" ]; then
  cat > activate_pipeline.sh << 'EOF'
#!/bin/bash
# Activate the Kelp GWAS pipeline environment

if command -v conda &> /dev/null; then
  eval "$(conda shell.bash hook)"
  conda activate kelp_gwas
  echo "Kelp GWAS environment activated!"
  echo "Python: $(which python)"
  echo "R: $(which R)"
else
  echo "Conda not found. Please install Miniconda or Anaconda."
  exit 1
fi
EOF

  chmod +x activate_pipeline.sh

  echo "Created activation script: activate_pipeline.sh"
  echo "Run: source activate_pipeline.sh"
  echo
fi

echo "============================================================================="
echo "Installation script completed!"
echo "============================================================================="
