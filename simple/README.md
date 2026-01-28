# Simple LOGAN-Style WGS (Local Snakemake)

This folder provides a **local, single-host WGS pipeline** in Snakemake that mirrors the core LOGAN flow while removing HPC/Slurm scatter/gather steps. It is intended for learning and small-scale runs on a single Ubuntu container.

## What this pipeline does (recommended single path)

1. **QC on raw reads**
   - `fastp` trimming + reports
   - `fastqc`, `fastq_screen`, `kraken2` (optional, but enabled by default)

2. **Alignment + post-processing**
   - `bwa-mem2` alignment
   - `samblaster` duplicate marking
   - `samtools` sort + index

3. **Base Quality Score Recalibration (BQSR)**
   - `gatk BaseRecalibrator` + `ApplyBQSR`

4. **Somatic SNV (tumor/normal)**
   - `gatk Mutect2` + contamination + read-orientation model + filtering

5. **Structural Variants (tumor/normal)**
   - `manta` somatic SV calling

6. **Copy Number Variants (tumor/normal)**
   - `cnvkit` batch mode (WGS)

7. **Reporting**
   - `multiqc` collects QC outputs

This is the **most recommended single path** for local use that tracks LOGAN's default best-practice choices without the full caller ensemble.

## Other LOGAN options (not enabled in this simplified path)

LOGAN can run many additional callers and filters. You can extend the Snakefile to add these if needed:

- **Somatic SNV**: VarDict 1.8.3, VarScan 2.4.4, Octopus (latest), Strelka 2.9.10, MuSE 2.x, LoFreq 0.0.1, DeepSomatic 1.8.0, Sage 3.4.4
- **SV**: GRIDSS 2.13.2 + GRIPSS 2.3.4, SvABA 1.2.0, SURVIVOR merge + AnnotSV 3.4.2
- **CNV**: FREEC 11.6, Sequenza utils 3.0.0, ASCAT (latest), CNVkit (git), HMF Purple 4.0.2 + Amber 4.0.1 + Cobalt 1.16
- **FFPE bias filtering**: SOBDetector (LOGAN ffpe module)
- **Annotation**: vcf2maf 102.0.0 + Ensembl VEP cache

These are listed here so you can expand the pipeline while keeping version parity with LOGAN.

## Files in this folder

- `Snakefile`: pipeline logic
- `config.json`: sample IDs
- `install.sh`: installs tools (LOGAN-matched versions)

## Setup

1. **Install tools** (inside your Ubuntu container):

```bash
bash install.sh
```

2. **Place input FASTQs** in this folder:
   - `{sample}_1.fq.gz` and `{sample}_2.fq.gz`
   - `{sample}` must match `TUMOR_SAMPLE` and `NORMAL_SAMPLE` in `config.json`

3. **Update `config.json`** with sample IDs:
   - `TUMOR_SAMPLE`, `NORMAL_SAMPLE`, `PAIR_ID`

`config.json` format (example):

```json
{
  "TUMOR_SAMPLE": "PTX901_T_wgs_lane1",
  "NORMAL_SAMPLE": "PTX901_N_blood_wgs_lane1",
  "PAIR_ID": "PTX901_T_vs_PTX901_N"
}
```

## Run (local)

```bash
snakemake --cores 8
```

Optional dry-run:

```bash
snakemake -n
```

## Optional features

- DeepVariant and annotation callers are not wired in this Snakefile. Add them if you want a broader LOGAN-style caller set.

## Output layout (default)

- `results/fastp/` - trimmed FASTQs + JSON/HTML
- `results/align/` - BAMs + BQSR BAMs
- `results/qc/` - fastqc, fastq_screen, kraken, mosdepth, qualimap, flagstat
- `results/snv/` - Mutect2 VCFs
- `results/sv/` - Manta somatic SVs
- `results/cnv/` - CNVkit outputs
- `results/multiqc/` - MultiQC report

## Version parity with LOGAN

These versions are matched to LOGAN container definitions:

- bwa-mem2 2.2.1
- samtools 1.20
- bcftools 1.20
- htslib 1.20
- samblaster 0.1.26
- fastp 0.24.0
- GATK 4.6.1.0
- Picard 3.2.0
- FastQC 0.12.1
- FastQ Screen 0.15.3
- Bowtie2 2.5.3
- Qualimap 2.3
- Mosdepth 0.3.8
- Kraken2 2.1.3
- KronaTools 2.8.1
- Manta 1.6.0
- MultiQC 1.23

CNVkit is installed from the upstream git repo, matching LOGAN's CNV container behavior (no pinned tag).

## Notes for learners

- LOGAN uses interval scatter/gather to scale on clusters. This local version runs **whole-genome steps** in one chunk for simplicity.
- The Mutect2 path here is the closest to LOGAN's tumor/normal best practice. It includes contamination modeling and read-orientation bias steps.
- Keep reference files consistent (same build for fasta, known-sites, germline resources, PON).
- `fastq_screen` needs a Bowtie2-based `fastq_screen.conf` with valid genome indexes.
- `kraken2` needs a pre-built database at `resources/hg38/kraken_db` (adjust the Snakefile if you store it elsewhere).
