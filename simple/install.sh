#!/usr/bin/env bash
set -euo pipefail

TOOLS_DIR="/opt/tools"
BIN_DIR="${TOOLS_DIR}/bin"
SRC_DIR="${TOOLS_DIR}/src"

mkdir -p "${BIN_DIR}" "${SRC_DIR}"

export PATH="${BIN_DIR}:${PATH}"

apt-get update
DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  bc \
  bzip2 \
  build-essential \
  ca-certificates \
  curl \
  git \
  libbz2-dev \
  libcurl4-openssl-dev \
  libgd-perl \
  liblzma-dev \
  libncurses5-dev \
  libncursesw5-dev \
  libssl-dev \
  perl \
  python3 \
  python3-pip \
  unzip \
  wget \
  zlib1g-dev \
  openjdk-17-jdk

# Optional: python2 for tools that still require it (Manta 1.6.0)
if ! command -v python2 >/dev/null 2>&1; then
  DEBIAN_FRONTEND=noninteractive apt-get install -y python2 || true
fi

# BWA-MEM2 2.2.1
cd "${SRC_DIR}"
wget -q https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2
rm -rf bwa-mem2-2.2.1_x64-linux
bunzip2 -c bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar -xf -
ln -sf "${SRC_DIR}/bwa-mem2-2.2.1_x64-linux/bwa-mem2" "${BIN_DIR}/bwa-mem2"
ln -sf "${SRC_DIR}/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2" "${BIN_DIR}/bwa-mem2.avx2"
ln -sf "${SRC_DIR}/bwa-mem2-2.2.1_x64-linux/bwa-mem2.sse41" "${BIN_DIR}/bwa-mem2.sse41"

# HTSlib/Samtools/Bcftools 1.20
cd "${SRC_DIR}"
wget -q https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
rm -rf htslib-1.20
tar -xjf htslib-1.20.tar.bz2
cd htslib-1.20
make -j"$(nproc)"
ln -sf "${SRC_DIR}/htslib-1.20" "${TOOLS_DIR}/htslib"

cd "${SRC_DIR}"
wget -q https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2
rm -rf samtools-1.20
tar -xjf samtools-1.20.tar.bz2
cd samtools-1.20
make -j"$(nproc)"
ln -sf "${SRC_DIR}/samtools-1.20/samtools" "${BIN_DIR}/samtools"

cd "${SRC_DIR}"
wget -q https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2
rm -rf bcftools-1.20
tar -xjf bcftools-1.20.tar.bz2
cd bcftools-1.20
make -j"$(nproc)"
ln -sf "${SRC_DIR}/bcftools-1.20/bcftools" "${BIN_DIR}/bcftools"

# Samblaster 0.1.26
cd "${SRC_DIR}"
rm -rf samblaster
git clone https://github.com/GregoryFaust/samblaster.git
cd samblaster
make -j"$(nproc)"
ln -sf "${SRC_DIR}/samblaster/samblaster" "${BIN_DIR}/samblaster"

# fastp 0.24.0
cd "${SRC_DIR}"
wget -q -O fastp https://github.com/OpenGene/fastp/releases/download/v0.24.0/fastp.0.24.0
chmod +x fastp
ln -sf "${SRC_DIR}/fastp" "${BIN_DIR}/fastp"

# GATK 4.6.1.0
cd "${SRC_DIR}"
wget -q https://github.com/broadinstitute/gatk/releases/download/4.6.1.0/gatk-4.6.1.0.zip
rm -rf gatk-4.6.1.0
unzip -q gatk-4.6.1.0.zip
ln -sf "${SRC_DIR}/gatk-4.6.1.0/gatk" "${BIN_DIR}/gatk"

# Picard 3.2.0
cd "${SRC_DIR}"
mkdir -p picard
wget -q -O picard/picard.jar https://github.com/broadinstitute/picard/releases/download/3.2.0/picard.jar

# FastQC 0.12.1
cd "${SRC_DIR}"
wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
rm -rf FastQC
unzip -q fastqc_v0.12.1.zip
chmod +x FastQC/fastqc
ln -sf "${SRC_DIR}/FastQC/fastqc" "${BIN_DIR}/fastqc"

# FastQ Screen 0.15.3
cd "${SRC_DIR}"
yes | perl -MCPAN -e "install GD"
yes | perl -MCPAN -e "install GD::Graph"
wget -q https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v0.15.3.tar.gz
rm -rf FastQ-Screen-0.15.3
tar -xzf v0.15.3.tar.gz
ln -sf "${SRC_DIR}/FastQ-Screen-0.15.3/fastq_screen" "${BIN_DIR}/fastq_screen"

# Bowtie2 2.5.3 (for FastQ Screen)
cd "${SRC_DIR}"
wget -q https://github.com/BenLangmead/bowtie2/releases/download/v2.5.3/bowtie2-2.5.3-linux-x86_64.zip
rm -rf bowtie2-2.5.3-linux-x86_64
unzip -q bowtie2-2.5.3-linux-x86_64.zip
ln -sf "${SRC_DIR}/bowtie2-2.5.3-linux-x86_64/bowtie2" "${BIN_DIR}/bowtie2"
ln -sf "${SRC_DIR}/bowtie2-2.5.3-linux-x86_64/bowtie2-build" "${BIN_DIR}/bowtie2-build"

# Qualimap 2.3
cd "${SRC_DIR}"
wget -q https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.3.zip
rm -rf qualimap_v2.3
unzip -q qualimap_v2.3.zip
ln -sf "${SRC_DIR}/qualimap_v2.3/qualimap" "${BIN_DIR}/qualimap"

# Mosdepth 0.3.8
cd "${SRC_DIR}"
wget -q -O mosdepth https://github.com/brentp/mosdepth/releases/download/v0.3.8/mosdepth
chmod +x mosdepth
ln -sf "${SRC_DIR}/mosdepth" "${BIN_DIR}/mosdepth"

# Kraken2 2.1.3
cd "${SRC_DIR}"
rm -rf kraken2-2.1.3
wget -q https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.3.tar.gz
tar -xzf v2.1.3.tar.gz
cd kraken2-2.1.3
./install_kraken2.sh "${SRC_DIR}/kraken2-2.1.3"
ln -sf "${SRC_DIR}/kraken2-2.1.3/kraken2" "${BIN_DIR}/kraken2"
ln -sf "${SRC_DIR}/kraken2-2.1.3/kraken2-build" "${BIN_DIR}/kraken2-build"

# KronaTools 2.8.1
cd "${SRC_DIR}"
rm -rf KronaTools-2.8.1
wget -q https://github.com/marbl/Krona/releases/download/v2.8.1/KronaTools-2.8.1.tar
tar -xf KronaTools-2.8.1.tar
cd KronaTools-2.8.1
./install.pl --prefix "${SRC_DIR}/KronaTools-2.8.1"
./updateTaxonomy.sh
ln -sf "${SRC_DIR}/KronaTools-2.8.1/bin/ktImportTaxonomy" "${BIN_DIR}/ktImportTaxonomy"

# Manta 1.6.0
cd "${SRC_DIR}"
rm -rf manta-1.6.0.centos6_x86_64
wget -q https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2
tar -xjf manta-1.6.0.centos6_x86_64.tar.bz2
ln -sf "${SRC_DIR}/manta-1.6.0.centos6_x86_64/bin/configManta.py" "${BIN_DIR}/configManta.py"

# CNVkit (LOGAN CNV container uses a git clone without a pinned version)
cd "${SRC_DIR}"
if [ ! -d "cnvkit" ]; then
  git clone https://github.com/etal/cnvkit
fi
cd cnvkit
pip3 install --upgrade pip
pip3 install -e .

# MultiQC 1.23
pip3 install multiqc==1.23

# Make PATH persistent for interactive shells
cat <<PROFILE > /etc/profile.d/logan_simple.sh
export PATH="${BIN_DIR}:\$PATH"
PROFILE
