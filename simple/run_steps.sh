#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Log major step completion to the results directory.
LOG_DIR="results"
LOG_FILE="${LOG_DIR}/run_steps.log"
mkdir -p "${LOG_DIR}"
touch "${LOG_FILE}"
log_step() {
  printf '%s\t%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "$1" >> "${LOG_FILE}"
}

# Sample IDs used to build FASTQ inputs and output prefixes.
CONFIG_JSON="${CONFIG_JSON:-${SCRIPT_DIR}/config.json}"
readarray -t SAMPLE_CONFIG < <(
  python3 - <<'PY' "${CONFIG_JSON}"
import json
import sys

with open(sys.argv[1]) as handle:
    config = json.load(handle)

print(config["TUMOR_SAMPLE"])
print(config["NORMAL_SAMPLE"])
print(config["PAIR_ID"])
PY
)
TUMOR_BASE_NAME="${SAMPLE_CONFIG[0]}"
NORMAL_BASE_NAME="${SAMPLE_CONFIG[1]}"
PAIR_ID="${SAMPLE_CONFIG[2]}"
# Thread count for tools that accept parallelism flags.
THREADS=8
TMPDIR="${PWD}/tmp"
mkdir -p "${TMPDIR}"

# Use AVX2-optimized bwa-mem2.
BWA_MEM2_BIN="bwa-mem2.avx2"

# Reference bundle layout expected under resources/hg38.
REF_DIR="${PWD}/resources/hg38"
REF_FASTA="${REF_DIR}/Homo_sapiens_assembly38.fasta"
REF_DICT="${REF_DIR}/Homo_sapiens_assembly38.dict"
KNOWN_SITES_DBSNP="${REF_DIR}/known_sites/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
KNOWN_SITES_MILLS="${REF_DIR}/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KNOWN_SITES_PHASE1="${REF_DIR}/known_sites/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
GERMLINE_RESOURCE="${REF_DIR}/somatic/af-only-gnomad.hg38.vcf.gz"
PON="${REF_DIR}/somatic/1000g_pon.hg38.vcf.gz"
INTERVALS="/home/ubuntu/LOGAN/assets/hg38_v0_wgs_calling_regions.hg38.interval_list"
CNVKIT_ACCESS="${REF_DIR}/cnvkit/access-10kb.hg38.bed"

# Optional contamination-screening resources.
FASTQ_SCREEN_CONF="${REF_DIR}/fastq_screen.conf"
KRAKEN_DB="${REF_DIR}/kraken_db"

# Upstream sources for reference files when missing.
REF_FASTA_URL="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
DBSNP_URL="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
DBSNP_TBI_URL="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"
MILLS_URL="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
MILLS_TBI_URL="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
PHASE1_URL="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
PHASE1_TBI_URL="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"
GERMLINE_URL="https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
GERMLINE_TBI_URL="https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi"
PON_URL="https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
PON_TBI_URL="https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi"

# Idempotent download helper for reference assets.
fetch() {
  local url="$1"
  local dest="$2"
  if [ ! -s "${dest}" ]; then
    mkdir -p "$(dirname "${dest}")"
    echo "Downloading $(basename "${dest}")..."
    curl -L -o "${dest}" "${url}"
  fi
}

# Ensure reference FASTA, indexes, and known-sites data are present.
prepare_resources() {
  mkdir -p "${REF_DIR}" "${REF_DIR}/known_sites" "${REF_DIR}/somatic" "${REF_DIR}/cnvkit"
  fetch "${REF_FASTA_URL}" "${REF_FASTA}"
  if [ ! -s "${REF_FASTA}.fai" ]; then
    samtools faidx "${REF_FASTA}"
  fi
  if [ ! -s "${REF_DICT}" ]; then
    gatk CreateSequenceDictionary -R "${REF_FASTA}" -O "${REF_DICT}"
  fi
  if [ ! -s "${REF_FASTA}.0123" ] || [ ! -s "${REF_FASTA}.bwt.2bit.64" ] || [ ! -s "${REF_FASTA}.sa" ]; then
    rm -f "${REF_FASTA}".{0123,amb,ann,bwt.2bit.64,pac,sa}
    "${BWA_MEM2_BIN}" index "${REF_FASTA}"
  fi
  fetch "${DBSNP_URL}" "${KNOWN_SITES_DBSNP}"
  fetch "${DBSNP_TBI_URL}" "${KNOWN_SITES_DBSNP}.tbi"
  fetch "${MILLS_URL}" "${KNOWN_SITES_MILLS}"
  fetch "${MILLS_TBI_URL}" "${KNOWN_SITES_MILLS}.tbi"
  fetch "${PHASE1_URL}" "${KNOWN_SITES_PHASE1}"
  fetch "${PHASE1_TBI_URL}" "${KNOWN_SITES_PHASE1}.tbi"
  fetch "${GERMLINE_URL}" "${GERMLINE_RESOURCE}"
  fetch "${GERMLINE_TBI_URL}" "${GERMLINE_RESOURCE}.tbi"
  fetch "${PON_URL}" "${PON}"
  fetch "${PON_TBI_URL}" "${PON}.tbi"
  if [ ! -s "${CNVKIT_ACCESS}" ]; then
    cnvkit.py access "${REF_FASTA}" -s 10000 -o "${CNVKIT_ACCESS}"
  fi
}

# Prepare reference data before any downstream processing.
prepare_resources
log_step "resources prepared"

# Step 1: Trim adapters/low-quality bases with fastp.
# Expects raw inputs named ${SAMPLE}_1.fq.gz / ${SAMPLE}_2.fq.gz in the current directory.
# Writes trimmed FASTQs plus per-sample JSON/HTML reports under results/fastp.
mkdir -p results/fastp
fastp -w "${THREADS}" \
    -i ${TUMOR_BASE_NAME}_1.fq.gz -I ${TUMOR_BASE_NAME}_2.fq.gz \
    -o results/fastp/${TUMOR_BASE_NAME}.R1.trimmed.fastq.gz -O results/fastp/${TUMOR_BASE_NAME}.R2.trimmed.fastq.gz \
    --json results/fastp/${TUMOR_BASE_NAME}.fastp.json \
    --html results/fastp/${TUMOR_BASE_NAME}.fastp.html
fastp -w "${THREADS}" \
    -i ${NORMAL_BASE_NAME}_1.fq.gz -I ${NORMAL_BASE_NAME}_2.fq.gz \
    -o results/fastp/${NORMAL_BASE_NAME}.R1.trimmed.fastq.gz -O results/fastp/${NORMAL_BASE_NAME}.R2.trimmed.fastq.gz \
    --json results/fastp/${NORMAL_BASE_NAME}.fastp.json \
    --html results/fastp/${NORMAL_BASE_NAME}.fastp.html
log_step "fastp complete"

emailme

# Purpose: Generate per-sample FastQC reports on trimmed reads.
mkdir -p results/qc/fastqc
fastqc -t "${THREADS}" -o results/qc/fastqc \
    results/fastp/${TUMOR_BASE_NAME}.R1.trimmed.fastq.gz results/fastp/${TUMOR_BASE_NAME}.R2.trimmed.fastq.gz
fastqc -t "${THREADS}" -o results/qc/fastqc \
    results/fastp/${NORMAL_BASE_NAME}.R1.trimmed.fastq.gz results/fastp/${NORMAL_BASE_NAME}.R2.trimmed.fastq.gz
log_step "fastqc complete"

emailme

# Purpose: Screen trimmed reads for contaminant genomes.
mkdir -p results/qc/fastq_screen
fastq_screen --conf "${FASTQ_SCREEN_CONF}" \
    --outdir results/qc/fastq_screen \
    --threads "${THREADS}" \
    --subset 1000000 \
    --aligner bowtie2 \
    --force \
    results/fastp/${TUMOR_BASE_NAME}.R1.trimmed.fastq.gz results/fastp/${TUMOR_BASE_NAME}.R2.trimmed.fastq.gz
fastq_screen --conf "${FASTQ_SCREEN_CONF}" \
    --outdir results/qc/fastq_screen \
    --threads "${THREADS}" \
    --subset 1000000 \
    --aligner bowtie2 \
    --force \
    results/fastp/${NORMAL_BASE_NAME}.R1.trimmed.fastq.gz results/fastp/${NORMAL_BASE_NAME}.R2.trimmed.fastq.gz
log_step "fastq_screen complete"

# Purpose: Classify reads taxonomically and build Krona summaries.
mkdir -p results/qc/kraken
kraken2 --db "${KRAKEN_DB}" \
    --threads "${THREADS}" --report results/qc/kraken/${TUMOR_BASE_NAME}.trimmed.kraken_bacteria.taxa.txt \
    --output - \
    --gzip-compressed \
    --paired results/fastp/${TUMOR_BASE_NAME}.R1.trimmed.fastq.gz results/fastp/${TUMOR_BASE_NAME}.R2.trimmed.fastq.gz
cut -f2,3 results/qc/kraken/${TUMOR_BASE_NAME}.trimmed.kraken_bacteria.taxa.txt | \
    ktImportTaxonomy - -o results/qc/kraken/${TUMOR_BASE_NAME}.trimmed.kraken_bacteria.krona.html
kraken2 --db "${KRAKEN_DB}" \
    --threads "${THREADS}" --report results/qc/kraken/${NORMAL_BASE_NAME}.trimmed.kraken_bacteria.taxa.txt \
    --output - \
    --gzip-compressed \
    --paired results/fastp/${NORMAL_BASE_NAME}.R1.trimmed.fastq.gz results/fastp/${NORMAL_BASE_NAME}.R2.trimmed.fastq.gz
cut -f2,3 results/qc/kraken/${NORMAL_BASE_NAME}.trimmed.kraken_bacteria.taxa.txt | \
    ktImportTaxonomy - -o results/qc/kraken/${NORMAL_BASE_NAME}.trimmed.kraken_bacteria.krona.html
log_step "kraken2 complete"

# Purpose: Align reads, mark duplicates, sort BAM, and write indexes.
mkdir -p results/align
mkdir -p "${TMPDIR}/bwa_${TUMOR_BASE_NAME}"
"${BWA_MEM2_BIN}" mem -M \
    -R '@RG\tID:${TUMOR_BASE_NAME}\tSM:${TUMOR_BASE_NAME}\tPL:illumina\tLB:${TUMOR_BASE_NAME}\tPU:${TUMOR_BASE_NAME}\tDS:wgs' \
    -t "${THREADS}" \
    "${REF_FASTA}" \
    results/fastp/${TUMOR_BASE_NAME}.R1.trimmed.fastq.gz results/fastp/${TUMOR_BASE_NAME}.R2.trimmed.fastq.gz | \
samblaster -M | \
samtools sort -T "${TMPDIR}/bwa_${TUMOR_BASE_NAME}/" \
    -@ "${THREADS}" -m 2G - \
    --write-index -o results/align/${TUMOR_BASE_NAME}.bam##idx##results/align/${TUMOR_BASE_NAME}.bam.bai

mkdir -p results/align
mkdir -p "${TMPDIR}/bwa_${NORMAL_BASE_NAME}"
"${BWA_MEM2_BIN}" mem -M \
    -R '@RG\tID:${NORMAL_BASE_NAME}\tSM:${NORMAL_BASE_NAME}\tPL:illumina\tLB:${NORMAL_BASE_NAME}\tPU:${NORMAL_BASE_NAME}\tDS:wgs' \
    -t "${THREADS}" \
    "${REF_FASTA}" \
    results/fastp/${NORMAL_BASE_NAME}.R1.trimmed.fastq.gz results/fastp/${NORMAL_BASE_NAME}.R2.trimmed.fastq.gz | \
samblaster -M | \
samtools sort -T "${TMPDIR}/bwa_${NORMAL_BASE_NAME}/" \
    -@ "${THREADS}" -m 2G - \
    --write-index -o results/align/${NORMAL_BASE_NAME}.bam##idx##results/align/${NORMAL_BASE_NAME}.bam.bai
log_step "alignment complete"

emailme

# Purpose: Build base quality recalibration tables.
mkdir -p results/bqsr
gatk --java-options '-Xmx10g' BaseRecalibrator \
    --input results/align/${TUMOR_BASE_NAME}.bam \
    --reference "${REF_FASTA}" \
    --known-sites "${KNOWN_SITES_DBSNP}" \
    --known-sites "${KNOWN_SITES_MILLS}" \
    --known-sites "${KNOWN_SITES_PHASE1}" \
    --intervals "${INTERVALS}" \
    --output results/bqsr/${TUMOR_BASE_NAME}.recal_data.grp
gatk --java-options '-Xmx10g' BaseRecalibrator \
    --input results/align/${NORMAL_BASE_NAME}.bam \
    --reference "${REF_FASTA}" \
    --known-sites "${KNOWN_SITES_DBSNP}" \
    --known-sites "${KNOWN_SITES_MILLS}" \
    --known-sites "${KNOWN_SITES_PHASE1}" \
    --intervals "${INTERVALS}" \
    --output results/bqsr/${NORMAL_BASE_NAME}.recal_data.grp
log_step "bqsr tables complete"

emailme

# Purpose: Apply BQSR and index recalibrated BAMs.
gatk --java-options '-Xmx20g' ApplyBQSR \
    --reference "${REF_FASTA}" \
    --input results/align/${TUMOR_BASE_NAME}.bam \
    --bqsr-recal-file results/bqsr/${TUMOR_BASE_NAME}.recal_data.grp \
    --output results/align/${TUMOR_BASE_NAME}.bqsr.bam \
    --use-jdk-inflater \
    --use-jdk-deflater
samtools index -@ "${THREADS}" results/align/${TUMOR_BASE_NAME}.bqsr.bam results/align/${TUMOR_BASE_NAME}.bqsr.bam.bai
gatk --java-options '-Xmx20g' ApplyBQSR \
    --reference "${REF_FASTA}" \
    --input results/align/${NORMAL_BASE_NAME}.bam \
    --bqsr-recal-file results/bqsr/${NORMAL_BASE_NAME}.recal_data.grp \
    --output results/align/${NORMAL_BASE_NAME}.bqsr.bam \
    --use-jdk-inflater \
    --use-jdk-deflater
samtools index -@ "${THREADS}" results/align/${NORMAL_BASE_NAME}.bqsr.bam results/align/${NORMAL_BASE_NAME}.bqsr.bam.bai
log_step "apply bqsr complete"

emailme

# Purpose: Summarize alignment metrics with flagstat.
mkdir -p results/qc/align
samtools flagstat results/align/${TUMOR_BASE_NAME}.bqsr.bam > results/qc/align/${TUMOR_BASE_NAME}.samtools_flagstat.txt
samtools flagstat results/align/${NORMAL_BASE_NAME}.bqsr.bam > results/qc/align/${NORMAL_BASE_NAME}.samtools_flagstat.txt
log_step "flagstat complete"

emailme

# Purpose: Estimate coverage distribution across the genome.
mkdir -p results/qc/align
mosdepth -n --fast-mode --by 500 results/qc/align/${TUMOR_BASE_NAME} results/align/${TUMOR_BASE_NAME}.bqsr.bam -t "${THREADS}"
mv results/qc/align/${TUMOR_BASE_NAME}.mosdepth.region.dist.txt results/qc/align/${TUMOR_BASE_NAME}.mosdepth.region.dist.txt
mv results/qc/align/${TUMOR_BASE_NAME}.mosdepth.summary.txt results/qc/align/${TUMOR_BASE_NAME}.mosdepth.summary.txt
mv results/qc/align/${TUMOR_BASE_NAME}.regions.bed.gz results/qc/align/${TUMOR_BASE_NAME}.regions.bed.gz
mv results/qc/align/${TUMOR_BASE_NAME}.regions.bed.gz.csi results/qc/align/${TUMOR_BASE_NAME}.regions.bed.gz.csi
mosdepth -n --fast-mode --by 500 results/qc/align/${NORMAL_BASE_NAME} results/align/${NORMAL_BASE_NAME}.bqsr.bam -t "${THREADS}"
mv results/qc/align/${NORMAL_BASE_NAME}.mosdepth.region.dist.txt results/qc/align/${NORMAL_BASE_NAME}.mosdepth.region.dist.txt
mv results/qc/align/${NORMAL_BASE_NAME}.mosdepth.summary.txt results/qc/align/${NORMAL_BASE_NAME}.mosdepth.summary.txt
mv results/qc/align/${NORMAL_BASE_NAME}.regions.bed.gz results/qc/align/${NORMAL_BASE_NAME}.regions.bed.gz
mv results/qc/align/${NORMAL_BASE_NAME}.regions.bed.gz.csi results/qc/align/${NORMAL_BASE_NAME}.regions.bed.gz.csi
log_step "mosdepth complete"

emailme

# Purpose: Generate Qualimap alignment QC reports.
mkdir -p results/qc/align
unset DISPLAY
qualimap bamqc -bam results/align/${TUMOR_BASE_NAME}.bqsr.bam \
    --java-mem-size=16G \
    -c -ip \
    -outdir results/qc/align/${TUMOR_BASE_NAME} \
    -outformat HTML \
    -nt "${THREADS}" \
    -nw 500 \
    -p NON-STRAND-SPECIFIC
mv results/qc/align/${TUMOR_BASE_NAME}/genome_results.txt results/qc/align/${TUMOR_BASE_NAME}_genome_results.txt
mv results/qc/align/${TUMOR_BASE_NAME}/qualimapReport.html results/qc/align/${TUMOR_BASE_NAME}_qualimapReport.html
qualimap bamqc -bam results/align/${NORMAL_BASE_NAME}.bqsr.bam \
    --java-mem-size=16G \
    -c -ip \
    -outdir results/qc/align/${NORMAL_BASE_NAME} \
    -outformat HTML \
    -nt "${THREADS}" \
    -nw 500 \
    -p NON-STRAND-SPECIFIC
mv results/qc/align/${NORMAL_BASE_NAME}/genome_results.txt results/qc/align/${NORMAL_BASE_NAME}_genome_results.txt
mv results/qc/align/${NORMAL_BASE_NAME}/qualimapReport.html results/qc/align/${NORMAL_BASE_NAME}_qualimapReport.html
log_step "qualimap complete"

emailme

# Purpose: Call somatic SNVs/indels with Mutect2.
mkdir -p results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2
gatk Mutect2 \
    --reference "${REF_FASTA}" \
    --intervals "${INTERVALS}" \
    --input results/align/${TUMOR_BASE_NAME}.bqsr.bam \
    --input results/align/${NORMAL_BASE_NAME}.bqsr.bam \
    --normal-sample ${NORMAL_BASE_NAME} \
    --tumor-sample ${TUMOR_BASE_NAME} \
    --germline-resource "${GERMLINE_RESOURCE}" \
    --panel-of-normals "${PON}" \
    --output results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.vcf.gz \
    --f1r2-tar-gz results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.f1r2.tar.gz \
    --independent-mates
log_step "mutect2 complete"

emailme

# Purpose: Collect pileups for contamination estimation.
mkdir -p results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2
gatk --java-options -Xmx48g GetPileupSummaries \
    -I results/align/${TUMOR_BASE_NAME}.bqsr.bam \
    -V "${GERMLINE_RESOURCE}" \
    --intervals "${INTERVALS}" \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.tpileup.table
gatk --java-options -Xmx48g GetPileupSummaries \
    -I results/align/${NORMAL_BASE_NAME}.bqsr.bam \
    -V "${GERMLINE_RESOURCE}" \
    --intervals "${INTERVALS}" \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.npileup.table
log_step "pileup summaries complete"

emailme

# Purpose: Aggregate pileups and compute contamination.
mkdir -p results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2
gatk GatherPileupSummaries \
    --sequence-dictionary "${REF_DICT}" \
    -I results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.tpileup.table \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.allpileups.table
gatk GatherPileupSummaries \
    --sequence-dictionary "${REF_DICT}" \
    -I results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.npileup.table \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.normal.allpileups.table
gatk CalculateContamination \
    -I results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.allpileups.table \
    --matched-normal results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.normal.allpileups.table \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.contamination.table
gatk CalculateContamination \
    -I results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.normal.allpileups.table \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.normal.contamination.table
log_step "contamination metrics complete"

emailme

# Purpose: Learn read-orientation bias model for filtering.
mkdir -p results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2
gatk LearnReadOrientationModel \
    --output results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.read-orientation-model.tar.gz \
    --input results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.f1r2.tar.gz
log_step "read orientation model complete"

emailme

# Purpose: Merge Mutect2 stats for downstream filtering.
mkdir -p results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2
gatk MergeMutectStats \
    --stats results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.vcf.gz.stats \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.final.stats
log_step "mutect2 stats merged"

emailme

# Purpose: Filter somatic calls and keep passing variants.
mkdir -p results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2
gatk SortVcf \
    -I results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.vcf.gz \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.concat.vcf.gz \
    --CREATE_INDEX
gatk FilterMutectCalls \
    -R "${REF_FASTA}" \
    -V results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.concat.vcf.gz \
    --ob-priors results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.read-orientation-model.tar.gz \
    --contamination-table results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.contamination.table \
    --stats results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.final.stats \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.marked.vcf.gz
gatk SelectVariants \
    -R "${REF_FASTA}" \
    --variant results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.marked.vcf.gz \
    --exclude-filtered \
    --output results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.final.vcf.gz
log_step "mutect2 filtering complete"

emailme

# Purpose: Normalize VCF representation and create an index.
bcftools sort results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.final.vcf.gz | \
    bcftools norm --threads "${THREADS}" --check-ref s \
    -f "${REF_FASTA}" -O v | \
    awk '{gsub(/\y[W|K|Y|R|S|M|B|D|H|V]\y/,"N",$4); OFS = "\t"; print}' | \
    sed '/^$/d' | \
    bcftools view - -Oz -o results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.norm.vcf.gz
bcftools index -t results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.norm.vcf.gz
log_step "mutect2 normalization complete"

emailme

# Purpose: Call somatic structural variants with Manta.
MANTA_CONFIG="$(command -v configManta.py)"
MANTA_CONFIG="$(readlink -f "${MANTA_CONFIG}")"
export PYTHONPATH="$(dirname "${MANTA_CONFIG}")${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p results/sv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/manta/run
"${MANTA_CONFIG}" \
    --tumorBam results/align/${TUMOR_BASE_NAME}.bqsr.bam \
    --normalBam results/align/${NORMAL_BASE_NAME}.bqsr.bam \
    --referenceFasta "${REF_FASTA}" \
    --runDir results/sv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/manta/run
results/sv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/manta/run/runWorkflow.py -m local -j "${THREADS}"
cp results/sv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/manta/run/results/variants/somaticSV.vcf.gz \
    results/sv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/manta/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.somaticSV.vcf.gz
cp results/sv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/manta/run/results/variants/somaticSV.vcf.gz.tbi \
    results/sv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/manta/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.somaticSV.vcf.gz.tbi
emailme
log_step "manta complete"

# Purpose: Call copy-number alterations with CNVkit.
mkdir -p results/cnv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/cnvkit
cnvkit.py batch \
    results/align/${TUMOR_BASE_NAME}.bqsr.bam -n results/align/${NORMAL_BASE_NAME}.bqsr.bam \
    -f "${REF_FASTA}" \
    -m wgs \
    --access "${CNVKIT_ACCESS}" \
    --output-reference results/cnv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/cnvkit/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.cnn \
    --output-dir results/cnv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/cnvkit
cp results/cnv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/cnvkit/${TUMOR_BASE_NAME}.bqsr.cns \
    results/cnv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/cnvkit/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.cns
log_step "cnvkit complete"

emailme

# Purpose: Aggregate QC metrics and reports across the run.
mkdir -p results/multiqc
multiqc results -o results/multiqc
log_step "multiqc complete"
