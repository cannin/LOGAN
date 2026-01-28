#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TMP_DIR="${SCRIPT_DIR}/tmp"

# Step 1: Choose a writable directory already listed on PATH for installs.
BIN_DIR=""
IFS=":" read -r -a PATH_DIRS <<< "${PATH}"
for path_dir in "${PATH_DIRS[@]}"; do
  if [ -z "${path_dir}" ]; then
    continue
  fi
  if [ -d "${path_dir}" ]; then
    if [ -w "${path_dir}" ]; then
      BIN_DIR="${path_dir}"
      break
    fi
  elif [ -w "$(dirname "${path_dir}")" ]; then
    mkdir -p "${path_dir}"
    BIN_DIR="${path_dir}"
    break
  fi
done

if [ -z "${BIN_DIR}" ]; then
  echo "No writable directory found in PATH; add one and rerun install.sh." >&2
  exit 1
fi

# Step 2: Use the repo-local temp directory for sources/build artifacts.
mkdir -p "${BIN_DIR}" "${TMP_DIR}"

# Step 3: Detect whether apt-get can be used (root or passwordless sudo).
APT_GET="apt-get"
if [ "$(id -u)" -ne 0 ] && command -v sudo >/dev/null 2>&1; then
  if sudo -n true 2>/dev/null; then
    APT_GET="sudo -n apt-get"
  else
    APT_GET=""
  fi
elif [ "$(id -u)" -ne 0 ]; then
  APT_GET=""
fi

# Step 4: Install OS-level dependencies when apt-get is available.
if [ -n "${APT_GET}" ]; then
  ${APT_GET} update
  DEBIAN_FRONTEND=noninteractive ${APT_GET} install -y --no-install-recommends \
  bc \
  bzip2 \
  build-essential \
  ca-certificates \
  curl \
  git \
  libbz2-dev \
  libcurl4-openssl-dev \
  libdbd-mysql-perl \
  libdbi-perl \
  libdeflate-dev \
  libgd-perl \
  libgd-graph-perl \
  libisal-dev \
  libffi-dev \
  liblist-moreutils-perl \
  libmodule-build-perl \
  libwww-perl \
  libsqlite3-dev \
  liblzma-dev \
  libncurses5-dev \
  libncursesw5-dev \
  libssl-dev \
  perl \
  python3 \
  python3-pip \
  python3-venv \
  unzip \
  wget \
  zlib1g-dev \
  openjdk-17-jdk
fi

# Step 5: Ensure python2 for legacy tools (build from source if absent).
# Optional: python2 for tools that still require it (Manta 1.6.0)
if [ -n "${APT_GET}" ] && ! command -v python2 >/dev/null 2>&1; then
  DEBIAN_FRONTEND=noninteractive ${APT_GET} install -y python2 || true
fi
if ! command -v python2 >/dev/null 2>&1; then
  PY2_DIR="${TMP_DIR}/python2"
  if [ ! -x "${PY2_DIR}/bin/python2" ]; then
    cd "${TMP_DIR}"
    rm -rf Python-2.7.18 Python-2.7.18.tgz
    wget -q -O Python-2.7.18.tgz https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz
    tar -xzf Python-2.7.18.tgz
    cd Python-2.7.18
    ./configure --prefix="${PY2_DIR}" LDFLAGS="-Wl,-rpath,${PY2_DIR}/lib"
    make -j"$(nproc)"
    make install
  fi
  ln -sf "${PY2_DIR}/bin/python2" "${BIN_DIR}/python2"
fi
# Provide a "python" entry if only python3 is available.
if ! command -v python >/dev/null 2>&1 && command -v python3 >/dev/null 2>&1; then
  ln -sf "$(command -v python3)" "${BIN_DIR}/python"
fi

# Step 6: Download/build bioinformatics tools and link binaries into BIN_DIR.
# BWA-MEM2 2.2.1
cd "${TMP_DIR}"
wget -q https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2
rm -rf bwa-mem2-2.2.1_x64-linux
bunzip2 -c bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar -xf -
ln -sf "${TMP_DIR}/bwa-mem2-2.2.1_x64-linux/bwa-mem2" "${BIN_DIR}/bwa-mem2"
ln -sf "${TMP_DIR}/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2" "${BIN_DIR}/bwa-mem2.avx2"
ln -sf "${TMP_DIR}/bwa-mem2-2.2.1_x64-linux/bwa-mem2.sse41" "${BIN_DIR}/bwa-mem2.sse41"

# HTSlib/Samtools/Bcftools 1.20
cd "${TMP_DIR}"
wget -q https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
rm -rf htslib-1.20
tar -xjf htslib-1.20.tar.bz2
cd htslib-1.20
make -j"$(nproc)"
ln -sf "${TMP_DIR}/htslib-1.20/tabix" "${BIN_DIR}/tabix"
ln -sf "${TMP_DIR}/htslib-1.20/bgzip" "${BIN_DIR}/bgzip"

cd "${TMP_DIR}"
wget -q https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2
rm -rf samtools-1.20
tar -xjf samtools-1.20.tar.bz2
cd samtools-1.20
make -j"$(nproc)"
ln -sf "${TMP_DIR}/samtools-1.20/samtools" "${BIN_DIR}/samtools"

cd "${TMP_DIR}"
wget -q https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2
rm -rf bcftools-1.20
tar -xjf bcftools-1.20.tar.bz2
cd bcftools-1.20
make -j"$(nproc)"
ln -sf "${TMP_DIR}/bcftools-1.20/bcftools" "${BIN_DIR}/bcftools"

# Samblaster 0.1.26
cd "${TMP_DIR}"
rm -rf samblaster
git clone https://github.com/GregoryFaust/samblaster.git
cd samblaster
make -j"$(nproc)"
ln -sf "${TMP_DIR}/samblaster/samblaster" "${BIN_DIR}/samblaster"

# fastp 0.24.0 (build from source; release has no binary assets)
cd "${TMP_DIR}"
rm -rf fastp-0.24.0 fastp-0.24.0.tar.gz fastp
wget -q -O fastp-0.24.0.tar.gz https://github.com/OpenGene/fastp/archive/refs/tags/v0.24.0.tar.gz
tar -xzf fastp-0.24.0.tar.gz
cd fastp-0.24.0
make -j"$(nproc)"
ln -sf "${TMP_DIR}/fastp-0.24.0/fastp" "${BIN_DIR}/fastp"

# GATK 4.6.1.0
cd "${TMP_DIR}"
rm -f gatk-4.6.1.0.zip
wget -q https://github.com/broadinstitute/gatk/releases/download/4.6.1.0/gatk-4.6.1.0.zip
rm -rf gatk-4.6.1.0
unzip -q gatk-4.6.1.0.zip
ln -sf "${TMP_DIR}/gatk-4.6.1.0/gatk" "${BIN_DIR}/gatk"

# Picard 3.2.0
cd "${TMP_DIR}"
mkdir -p picard
wget -q -O picard/picard.jar https://github.com/broadinstitute/picard/releases/download/3.2.0/picard.jar

# FastQC 0.12.1
cd "${TMP_DIR}"
wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
rm -rf FastQC
unzip -q fastqc_v0.12.1.zip
chmod +x FastQC/fastqc
ln -sf "${TMP_DIR}/FastQC/fastqc" "${BIN_DIR}/fastqc"

# FastQ Screen 0.15.3
cd "${TMP_DIR}"
if ! perl -MGD -e 1 >/dev/null 2>&1; then
  yes | PERL_MM_USE_DEFAULT=1 perl -MCPAN -e "install GD" || true
fi
if ! perl -MGD::Graph -e 1 >/dev/null 2>&1; then
  yes | PERL_MM_USE_DEFAULT=1 perl -MCPAN -e "install GD::Graph" || true
fi
wget -q https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v0.15.3.tar.gz
rm -rf FastQ-Screen-0.15.3
tar -xzf v0.15.3.tar.gz
ln -sf "${TMP_DIR}/FastQ-Screen-0.15.3/fastq_screen" "${BIN_DIR}/fastq_screen"

# Bowtie2 2.5.3 (for FastQ Screen)
cd "${TMP_DIR}"
wget -q https://github.com/BenLangmead/bowtie2/releases/download/v2.5.3/bowtie2-2.5.3-linux-x86_64.zip
rm -rf bowtie2-2.5.3-linux-x86_64
unzip -q bowtie2-2.5.3-linux-x86_64.zip
ln -sf "${TMP_DIR}/bowtie2-2.5.3-linux-x86_64/bowtie2" "${BIN_DIR}/bowtie2"
ln -sf "${TMP_DIR}/bowtie2-2.5.3-linux-x86_64/bowtie2-build" "${BIN_DIR}/bowtie2-build"

# Qualimap 2.3
cd "${TMP_DIR}"
curl -L -o qualimap_v2.3.zip https://api.bitbucket.org/2.0/repositories/kokonech/qualimap/downloads/qualimap_v2.3.zip
rm -rf qualimap_v2.3
unzip -q qualimap_v2.3.zip
ln -sf "${TMP_DIR}/qualimap_v2.3/qualimap" "${BIN_DIR}/qualimap"

# Mosdepth 0.3.8
cd "${TMP_DIR}"
wget -q -O mosdepth https://github.com/brentp/mosdepth/releases/download/v0.3.8/mosdepth
chmod +x mosdepth
ln -sf "${TMP_DIR}/mosdepth" "${BIN_DIR}/mosdepth"

# Kraken2 2.1.3
cd "${TMP_DIR}"
rm -rf kraken2-2.1.3
wget -q https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.3.tar.gz
tar -xzf v2.1.3.tar.gz
cd kraken2-2.1.3
./install_kraken2.sh "${TMP_DIR}/kraken2-2.1.3"
ln -sf "${TMP_DIR}/kraken2-2.1.3/kraken2" "${BIN_DIR}/kraken2"
ln -sf "${TMP_DIR}/kraken2-2.1.3/kraken2-build" "${BIN_DIR}/kraken2-build"

# KronaTools 2.8.1
cd "${TMP_DIR}"
rm -rf KronaTools-2.8.1
wget -q https://github.com/marbl/Krona/releases/download/v2.8.1/KronaTools-2.8.1.tar
tar -xf KronaTools-2.8.1.tar
cd KronaTools-2.8.1
./install.pl --prefix "${TMP_DIR}/KronaTools-2.8.1"
./updateTaxonomy.sh
ln -sf "${TMP_DIR}/KronaTools-2.8.1/bin/ktImportTaxonomy" "${BIN_DIR}/ktImportTaxonomy"

# Manta 1.6.0
cd "${TMP_DIR}"
rm -rf manta-1.6.0.centos6_x86_64
wget -q https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2
tar -xjf manta-1.6.0.centos6_x86_64.tar.bz2
ln -sf "${TMP_DIR}/manta-1.6.0.centos6_x86_64/bin/configManta.py" "${BIN_DIR}/configManta.py"

# vcf2maf (annotation helper; requires VEP cache/config at runtime)
cd "${TMP_DIR}"
if [ ! -d "vcf2maf" ]; then
  git clone https://github.com/mskcc/vcf2maf.git
else
  (cd vcf2maf && git fetch --all --prune)
fi
ln -sf "${TMP_DIR}/vcf2maf/vcf2maf.pl" "${BIN_DIR}/vcf2maf.pl"
ln -sf "${TMP_DIR}/vcf2maf/vcf2maf.pl" "${BIN_DIR}/vcf2maf"

# VEP (needed for vcf2maf; installs local cache under tmp/vep_cache)
VEP_TAG="${VEP_TAG:-release/115.0}"
VEP_CACHE_VERSION="${VEP_CACHE_VERSION:-115}"
VEP_DIR="${TMP_DIR}/vep"
VEP_CACHE_DIR="${TMP_DIR}/vep_cache"
VEP_API_DIR="${TMP_DIR}/vep_api"
mkdir -p "${VEP_API_DIR}"
if [ ! -x "${VEP_DIR}/vep" ] && [ ! -x "${VEP_DIR}/variant_effect_predictor.pl" ]; then
  cd "${TMP_DIR}"
  VEP_TARBALL="ensembl-vep-${VEP_TAG//\//-}.tar.gz"
  if [ ! -s "${VEP_TARBALL}" ]; then
    rm -f "${VEP_TARBALL}"
    curl -L -o "${VEP_TARBALL}" "https://github.com/Ensembl/ensembl-vep/archive/refs/tags/${VEP_TAG}.tar.gz"
  fi
  if [ ! -s "${VEP_TARBALL}" ]; then
    echo "Failed to download VEP tarball for tag ${VEP_TAG}." >&2
    exit 1
  fi
  rm -rf "${VEP_DIR}"
  mkdir -p "${VEP_DIR}"
  tar -xzf "${VEP_TARBALL}" -C "${VEP_DIR}" --strip-components=1
  cd "${VEP_DIR}"
  perl INSTALL.pl \
    --AUTO a \
    --NO_UPDATE \
    --NO_TEST \
    --DESTDIR "${VEP_API_DIR}"
fi
if [ ! -f "${VEP_API_DIR}/Bio/EnsEMBL/Registry.pm" ]; then
  cd "${VEP_DIR}"
  perl INSTALL.pl \
    --AUTO a \
    --NO_UPDATE \
    --NO_TEST \
    --DESTDIR "${VEP_API_DIR}"
fi
if [ ! -d "${VEP_CACHE_DIR}/homo_sapiens" ]; then
  mkdir -p "${VEP_CACHE_DIR}"
  cd "${VEP_DIR}"
  perl INSTALL.pl \
    --AUTO c \
    --NO_UPDATE \
    --DESTDIR "${VEP_API_DIR}" \
    --CACHEDIR "${VEP_CACHE_DIR}" \
    --SPECIES homo_sapiens \
    --ASSEMBLY GRCh38 \
    --CACHE_VERSION "${VEP_CACHE_VERSION}"
fi
if [ -x "${VEP_DIR}/vep" ]; then
  ln -sf "${VEP_DIR}/vep" "${BIN_DIR}/vep"
fi
if [ -x "${VEP_DIR}/variant_effect_predictor.pl" ]; then
  ln -sf "${VEP_DIR}/variant_effect_predictor.pl" "${BIN_DIR}/variant_effect_predictor.pl"
fi

# Step 7: Create a python venv and install python-based tools.
# CNVkit (LOGAN CNV container uses a git clone without a pinned version)
cd "${TMP_DIR}"
if [ ! -d "cnvkit" ]; then
  git clone https://github.com/etal/cnvkit
fi
VENV_DIR="${TMP_DIR}/venv"
if [ ! -d "${VENV_DIR}" ]; then
  python3 -m venv "${VENV_DIR}"
fi
cd cnvkit
"${VENV_DIR}/bin/pip" install --upgrade pip
"${VENV_DIR}/bin/pip" install -e .
ln -sf "${VENV_DIR}/bin/cnvkit.py" "${BIN_DIR}/cnvkit.py"

# MultiQC 1.23
"${VENV_DIR}/bin/pip" install multiqc==1.23
ln -sf "${VENV_DIR}/bin/multiqc" "${BIN_DIR}/multiqc"

# Snakemake (workflow runner)
"${VENV_DIR}/bin/pip" install snakemake==7.32.4
# PuLP 2.7.0 is required for Snakemake; newer releases removed list_solvers.
"${VENV_DIR}/bin/pip" install pulp==2.7.0
ln -sf "${VENV_DIR}/bin/snakemake" "${BIN_DIR}/snakemake"
