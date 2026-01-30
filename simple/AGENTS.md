# Simple Pipeline Guidelines

## Scope
This folder contains a local, single-host WGS pipeline implemented in `Snakefile`. It mirrors the core LOGAN flow without cluster scatter/gather. Use this guide when editing or running the simple pipeline only.

## Project Layout
- `Snakefile`: Snakemake workflow (all configuration is embedded here; no config file).
- `run_steps.sh`: Reference shell workflow (do not modify; keep changes in `Snakefile`).
- `install.sh`: Installs required tools into the local environment.
- `resources/hg38/`: Reference assets (FASTA, indexes, known-sites, optional QC databases).
- `results/`: Pipeline outputs and logs.

## Run Commands
- Install tools:
  - `bash install.sh`
- Dry run:
  - `snakemake -n --cores 8`
- Full run (example):
  - `snakemake --cores 30 --rerun-incomplete`
- Unlock after interruption:
  - `snakemake --unlock`

## Inputs and Outputs
- Inputs are raw FASTQs in the working directory:
  - `{sample}_1.fq.gz` and `{sample}_2.fq.gz`
- Sample IDs are defined in `Snakefile` (`TUMOR_SAMPLE`, `NORMAL_SAMPLE`, `PAIR_ID`). Update them there.
- Outputs are written under `results/` with subfolders for `fastp`, `qc`, `align`, `snv`, `sv`, `cnv`, and `multiqc`.
- Each rule appends a completion record to `results/rule_finish.log` and then calls `notifyme`.

## Tooling and References
- Reference layout follows `resources/hg38/` (FASTA, `.fai`, `.dict`, and bwa-mem2 indexes).
- `fastq_screen` reads `resources/hg38/fastq_screen.conf` and requires a Bowtie2 index.
- `kraken2` expects a database in `resources/hg38/kraken_db`.

## Resource Provenance (How They Are Retrieved)
- Resource paths are declared near the top of `Snakefile` (`REF_DIR`, `KNOWN_SITES_*`, `GERMLINE_RESOURCE`, `PON`, `FASTQ_SCREEN_CONF`, `KRAKEN_DB`, `CNVKIT_ACCESS`).
- Reference FASTA and known-sites VCFs mirror `run_steps.sh` sources (Broad GCP public references and GATK best practices). See the URLs in `simple/run_steps.sh` for exact downloads.
- bwa-mem2 indexes are generated locally: `bwa-mem2.avx2 index resources/hg38/Homo_sapiens_assembly38.fasta`.
- `CNVKIT_ACCESS` is created locally: `cnvkit.py access <fasta> -s 10000 -o resources/hg38/cnvkit/access-10kb.hg38.bed`.
- `fastq_screen` uses a Bowtie2 index built from the FASTA and referenced by `resources/hg38/fastq_screen.conf` (prefix `resources/hg38/hg38_bowtie2`).
- `kraken2` DB is built locally with `kraken2-build --download-taxonomy` + `--download-library bacteria` + `--build` under `resources/hg38/kraken_db` (use `--use-ftp` if rsync is blocked).

## Editing Conventions
- Keep changes simple and aligned with `run_steps.sh` behavior.
- Avoid adding new conditionals or complex logic in `Snakefile`.
- Use ASCII and keep shell blocks readable.
