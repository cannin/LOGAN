#!/usr/bin/env bash
set -euo pipefail

TUMOR_BASE_NAME="TUMOR_1"
NORMAL_BASE_NAME="NORMAL_1"

# Purpose: Trim and quality-filter raw reads.
module load fastp
mkdir -p results/fastp
fastp -w 8 \
    -i ${TUMOR_BASE_NAME}_1.fq.gz -I ${TUMOR_BASE_NAME}_2.fq.gz \
    -o results/fastp/${TUMOR_BASE_NAME}.R1.trimmed.fastq.gz -O results/fastp/${TUMOR_BASE_NAME}.R2.trimmed.fastq.gz \
    --json results/fastp/${TUMOR_BASE_NAME}.fastp.json \
    --html results/fastp/${TUMOR_BASE_NAME}.fastp.html
fastp -w 8 \
    -i ${NORMAL_BASE_NAME}_1.fq.gz -I ${NORMAL_BASE_NAME}_2.fq.gz \
    -o results/fastp/${NORMAL_BASE_NAME}.R1.trimmed.fastq.gz -O results/fastp/${NORMAL_BASE_NAME}.R2.trimmed.fastq.gz \
    --json results/fastp/${NORMAL_BASE_NAME}.fastp.json \
    --html results/fastp/${NORMAL_BASE_NAME}.fastp.html

# Purpose: Generate per-sample FastQC reports on trimmed reads.
module load fastqc
mkdir -p results/qc/fastqc
fastqc -t 8 -o results/qc/fastqc \
    results/fastp/${TUMOR_BASE_NAME}.R1.trimmed.fastq.gz results/fastp/${TUMOR_BASE_NAME}.R2.trimmed.fastq.gz
fastqc -t 8 -o results/qc/fastqc \
    results/fastp/${NORMAL_BASE_NAME}.R1.trimmed.fastq.gz results/fastp/${NORMAL_BASE_NAME}.R2.trimmed.fastq.gz

# Purpose: Screen trimmed reads for contaminant genomes.
module load fastq_screen bowtie2
mkdir -p results/qc/fastq_screen
fastq_screen --conf ../conf/fastq_screen.conf \
    --outdir results/qc/fastq_screen \
    --threads 8 \
    --subset 1000000 \
    --aligner bowtie2 \
    --force \
    results/fastp/${TUMOR_BASE_NAME}.R1.trimmed.fastq.gz results/fastp/${TUMOR_BASE_NAME}.R2.trimmed.fastq.gz
fastq_screen --conf ../conf/fastq_screen.conf \
    --outdir results/qc/fastq_screen \
    --threads 8 \
    --subset 1000000 \
    --aligner bowtie2 \
    --force \
    results/fastp/${NORMAL_BASE_NAME}.R1.trimmed.fastq.gz results/fastp/${NORMAL_BASE_NAME}.R2.trimmed.fastq.gz

# Purpose: Classify reads taxonomically and build Krona summaries.
module load kraken kronatools
mkdir -p results/qc/kraken
kraken2 --db /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/kraken/20180907_standard_kraken2 \
    --threads 8 --report results/qc/kraken/${TUMOR_BASE_NAME}.trimmed.kraken_bacteria.taxa.txt \
    --output - \
    --gzip-compressed \
    --paired results/fastp/${TUMOR_BASE_NAME}.R1.trimmed.fastq.gz results/fastp/${TUMOR_BASE_NAME}.R2.trimmed.fastq.gz
cut -f2,3 results/qc/kraken/${TUMOR_BASE_NAME}.trimmed.kraken_bacteria.taxa.txt | \
    ktImportTaxonomy - -o results/qc/kraken/${TUMOR_BASE_NAME}.trimmed.kraken_bacteria.krona.html
kraken2 --db /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/kraken/20180907_standard_kraken2 \
    --threads 8 --report results/qc/kraken/${NORMAL_BASE_NAME}.trimmed.kraken_bacteria.taxa.txt \
    --output - \
    --gzip-compressed \
    --paired results/fastp/${NORMAL_BASE_NAME}.R1.trimmed.fastq.gz results/fastp/${NORMAL_BASE_NAME}.R2.trimmed.fastq.gz
cut -f2,3 results/qc/kraken/${NORMAL_BASE_NAME}.trimmed.kraken_bacteria.taxa.txt | \
    ktImportTaxonomy - -o results/qc/kraken/${NORMAL_BASE_NAME}.trimmed.kraken_bacteria.krona.html

# Purpose: Align reads, mark duplicates, sort BAM, and write indexes.
module load bwa-mem2 samblaster samtools
mkdir -p results/align
mkdir -p /tmp/logan_simple/bwa_${TUMOR_BASE_NAME}
bwa-mem2.avx2 mem -M \
    -R '@RG\tID:${TUMOR_BASE_NAME}\tSM:${TUMOR_BASE_NAME}\tPL:illumina\tLB:${TUMOR_BASE_NAME}\tPU:${TUMOR_BASE_NAME}\tDS:wgs' \
    -t 8 \
    /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    results/fastp/${TUMOR_BASE_NAME}.R1.trimmed.fastq.gz results/fastp/${TUMOR_BASE_NAME}.R2.trimmed.fastq.gz | \
samblaster -M | \
samtools sort -T /tmp/logan_simple/bwa_${TUMOR_BASE_NAME}/ \
    -@ 8 -m 2G - \
    --write-index -o results/align/${TUMOR_BASE_NAME}.bam##idx##results/align/${TUMOR_BASE_NAME}.bam.bai

mkdir -p results/align
mkdir -p /tmp/logan_simple/bwa_${NORMAL_BASE_NAME}
bwa-mem2.avx2 mem -M \
    -R '@RG\tID:${NORMAL_BASE_NAME}\tSM:${NORMAL_BASE_NAME}\tPL:illumina\tLB:${NORMAL_BASE_NAME}\tPU:${NORMAL_BASE_NAME}\tDS:wgs' \
    -t 8 \
    /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    results/fastp/${NORMAL_BASE_NAME}.R1.trimmed.fastq.gz results/fastp/${NORMAL_BASE_NAME}.R2.trimmed.fastq.gz | \
samblaster -M | \
samtools sort -T /tmp/logan_simple/bwa_${NORMAL_BASE_NAME}/ \
    -@ 8 -m 2G - \
    --write-index -o results/align/${NORMAL_BASE_NAME}.bam##idx##results/align/${NORMAL_BASE_NAME}.bam.bai

# Purpose: Build base quality recalibration tables.
module load GATK
mkdir -p results/bqsr
gatk --java-options '-Xmx10g' BaseRecalibrator \
    --input results/align/${TUMOR_BASE_NAME}.bam \
    --reference /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    --known-sites /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/GATK_resource_bundle/dbsnp_138.hg38.vcf.gz \
    --known-sites /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-sites /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/GATK_resource_bundle/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz \
    --intervals /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list \
    --output results/bqsr/${TUMOR_BASE_NAME}.recal_data.grp
gatk --java-options '-Xmx10g' BaseRecalibrator \
    --input results/align/${NORMAL_BASE_NAME}.bam \
    --reference /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    --known-sites /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/GATK_resource_bundle/dbsnp_138.hg38.vcf.gz \
    --known-sites /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-sites /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/GATK_resource_bundle/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz \
    --intervals /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list \
    --output results/bqsr/${NORMAL_BASE_NAME}.recal_data.grp

# Purpose: Apply BQSR and index recalibrated BAMs.
module load GATK samtools
gatk --java-options '-Xmx20g' ApplyBQSR \
    --reference /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    --input results/align/${TUMOR_BASE_NAME}.bam \
    --bqsr-recal-file results/bqsr/${TUMOR_BASE_NAME}.recal_data.grp \
    --output results/align/${TUMOR_BASE_NAME}.bqsr.bam \
    --use-jdk-inflater \
    --use-jdk-deflater
samtools index -@ 8 results/align/${TUMOR_BASE_NAME}.bqsr.bam results/align/${TUMOR_BASE_NAME}.bqsr.bam.bai
gatk --java-options '-Xmx20g' ApplyBQSR \
    --reference /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    --input results/align/${NORMAL_BASE_NAME}.bam \
    --bqsr-recal-file results/bqsr/${NORMAL_BASE_NAME}.recal_data.grp \
    --output results/align/${NORMAL_BASE_NAME}.bqsr.bam \
    --use-jdk-inflater \
    --use-jdk-deflater
samtools index -@ 8 results/align/${NORMAL_BASE_NAME}.bqsr.bam results/align/${NORMAL_BASE_NAME}.bqsr.bam.bai

# Purpose: Summarize alignment metrics with flagstat.
module load samtools
mkdir -p results/qc/align
samtools flagstat results/align/${TUMOR_BASE_NAME}.bqsr.bam > results/qc/align/${TUMOR_BASE_NAME}.samtools_flagstat.txt
samtools flagstat results/align/${NORMAL_BASE_NAME}.bqsr.bam > results/qc/align/${NORMAL_BASE_NAME}.samtools_flagstat.txt

# Purpose: Estimate coverage distribution across the genome.
module load mosdepth
mkdir -p results/qc/align
mosdepth -n --fast-mode --by 500 results/qc/align/${TUMOR_BASE_NAME} results/align/${TUMOR_BASE_NAME}.bqsr.bam -t 8
mv results/qc/align/${TUMOR_BASE_NAME}.mosdepth.region.dist.txt results/qc/align/${TUMOR_BASE_NAME}.mosdepth.region.dist.txt
mv results/qc/align/${TUMOR_BASE_NAME}.mosdepth.summary.txt results/qc/align/${TUMOR_BASE_NAME}.mosdepth.summary.txt
mv results/qc/align/${TUMOR_BASE_NAME}.regions.bed.gz results/qc/align/${TUMOR_BASE_NAME}.regions.bed.gz
mv results/qc/align/${TUMOR_BASE_NAME}.regions.bed.gz.csi results/qc/align/${TUMOR_BASE_NAME}.regions.bed.gz.csi
mosdepth -n --fast-mode --by 500 results/qc/align/${NORMAL_BASE_NAME} results/align/${NORMAL_BASE_NAME}.bqsr.bam -t 8
mv results/qc/align/${NORMAL_BASE_NAME}.mosdepth.region.dist.txt results/qc/align/${NORMAL_BASE_NAME}.mosdepth.region.dist.txt
mv results/qc/align/${NORMAL_BASE_NAME}.mosdepth.summary.txt results/qc/align/${NORMAL_BASE_NAME}.mosdepth.summary.txt
mv results/qc/align/${NORMAL_BASE_NAME}.regions.bed.gz results/qc/align/${NORMAL_BASE_NAME}.regions.bed.gz
mv results/qc/align/${NORMAL_BASE_NAME}.regions.bed.gz.csi results/qc/align/${NORMAL_BASE_NAME}.regions.bed.gz.csi

# Purpose: Generate Qualimap alignment QC reports.
module load qualimap
mkdir -p results/qc/align
unset DISPLAY
qualimap bamqc -bam results/align/${TUMOR_BASE_NAME}.bqsr.bam \
    --java-mem-size=96G \
    -c -ip \
    -outdir results/qc/align/${TUMOR_BASE_NAME} \
    -outformat HTML \
    -nt 8 \
    -nw 500 \
    -p NON-STRAND-SPECIFIC
mv results/qc/align/${TUMOR_BASE_NAME}/genome_results.txt results/qc/align/${TUMOR_BASE_NAME}_genome_results.txt
mv results/qc/align/${TUMOR_BASE_NAME}/qualimapReport.html results/qc/align/${TUMOR_BASE_NAME}_qualimapReport.html
qualimap bamqc -bam results/align/${NORMAL_BASE_NAME}.bqsr.bam \
    --java-mem-size=96G \
    -c -ip \
    -outdir results/qc/align/${NORMAL_BASE_NAME} \
    -outformat HTML \
    -nt 8 \
    -nw 500 \
    -p NON-STRAND-SPECIFIC
mv results/qc/align/${NORMAL_BASE_NAME}/genome_results.txt results/qc/align/${NORMAL_BASE_NAME}_genome_results.txt
mv results/qc/align/${NORMAL_BASE_NAME}/qualimapReport.html results/qc/align/${NORMAL_BASE_NAME}_qualimapReport.html

# Purpose: Call somatic SNVs/indels with Mutect2.
module load GATK
mkdir -p results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2
gatk Mutect2 \
    --reference /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    --intervals /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list \
    --input results/align/${TUMOR_BASE_NAME}.bqsr.bam \
    --input results/align/${NORMAL_BASE_NAME}.bqsr.bam \
    --normal-sample ${NORMAL_BASE_NAME} \
    --tumor-sample ${TUMOR_BASE_NAME} \
    --germline-resource /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/GATK_GRCh38.d1.vd1/somatic-hg38-af-only-gnomad.hg38.vcf.gz \
    --panel-of-normals /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/PON/MuTect2.PON.5210.vcf.gz \
    --output results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.vcf.gz \
    --f1r2-tar-gz results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.f1r2.tar.gz \
    --independent-mates

# Purpose: Collect pileups for contamination estimation.
module load GATK
mkdir -p results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2
gatk --java-options -Xmx48g GetPileupSummaries \
    -I results/align/${TUMOR_BASE_NAME}.bqsr.bam \
    -V /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/GATK_GRCh38.d1.vd1/somatic-hg38-af-only-gnomad.hg38.vcf.gz \
    --intervals /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.tpileup.table
gatk --java-options -Xmx48g GetPileupSummaries \
    -I results/align/${NORMAL_BASE_NAME}.bqsr.bam \
    -V /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/GATK_GRCh38.d1.vd1/somatic-hg38-af-only-gnomad.hg38.vcf.gz \
    --intervals /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.npileup.table

# Purpose: Aggregate pileups and compute contamination.
module load GATK
mkdir -p results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2
gatk GatherPileupSummaries \
    --sequence-dictionary /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem/GRCh38.d1.vd1.dict \
    -I results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.tpileup.table \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.allpileups.table
gatk GatherPileupSummaries \
    --sequence-dictionary /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem/GRCh38.d1.vd1.dict \
    -I results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.npileup.table \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.normal.allpileups.table
gatk CalculateContamination \
    -I results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.allpileups.table \
    --matched-normal results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.normal.allpileups.table \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.contamination.table
gatk CalculateContamination \
    -I results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.normal.allpileups.table \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.normal.contamination.table

# Purpose: Learn read-orientation bias model for filtering.
module load GATK
mkdir -p results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2
gatk LearnReadOrientationModel \
    --output results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.read-orientation-model.tar.gz \
    --input results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.f1r2.tar.gz

# Purpose: Merge Mutect2 stats for downstream filtering.
module load GATK
mkdir -p results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2
gatk MergeMutectStats \
    --stats results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.vcf.gz.stats \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.final.stats

# Purpose: Filter somatic calls and keep passing variants.
module load GATK
mkdir -p results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2
gatk SortVcf \
    -I results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.vcf.gz \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.concat.vcf.gz \
    --CREATE_INDEX
gatk FilterMutectCalls \
    -R /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    -V results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.concat.vcf.gz \
    --ob-priors results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.read-orientation-model.tar.gz \
    --contamination-table results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.contamination.table \
    --stats results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.final.stats \
    -O results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.marked.vcf.gz
gatk SelectVariants \
    -R /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    --variant results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.marked.vcf.gz \
    --exclude-filtered \
    --output results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.final.vcf.gz

# Purpose: Normalize VCF representation and create an index.
module load bcftools
bcftools sort results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.final.vcf.gz | \
    bcftools norm --threads 8 --check-ref s \
    -f /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa -O v | \
    awk '{gsub(/\y[W|K|Y|R|S|M|B|D|H|V]\y/,"N",$4); OFS = "\t"; print}' | \
    sed '/^$/d' | \
    bcftools view - -Oz -o results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.norm.vcf.gz
bcftools index -t results/snv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/mutect2/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.mutect2.norm.vcf.gz

# Purpose: Call somatic structural variants with Manta.
module load manta
mkdir -p results/sv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/manta/run
configManta.py \
    --tumorBam results/align/${TUMOR_BASE_NAME}.bqsr.bam \
    --normalBam results/align/${NORMAL_BASE_NAME}.bqsr.bam \
    --referenceFasta /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    --runDir results/sv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/manta/run
results/sv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/manta/run/runWorkflow.py -m local -j 8
cp results/sv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/manta/run/results/variants/somaticSV.vcf.gz \
    results/sv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/manta/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.somaticSV.vcf.gz
cp results/sv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/manta/run/results/variants/somaticSV.vcf.gz.tbi \
    results/sv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/manta/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.somaticSV.vcf.gz.tbi

# Purpose: Call copy-number alterations with CNVkit.
module load cnvkit
mkdir -p results/cnv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/cnvkit
cnvkit.py batch \
    results/align/${TUMOR_BASE_NAME}.bqsr.bam -n results/align/${NORMAL_BASE_NAME}.bqsr.bam \
    -f /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    -m wgs \
    --access /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/cnvkit/access-10kb.hg38.bed \
    --output-reference results/cnv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/cnvkit/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.cnn \
    --output-dir results/cnv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/cnvkit
cp results/cnv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/cnvkit/${TUMOR_BASE_NAME}.bqsr.cns \
    results/cnv/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}/cnvkit/${TUMOR_BASE_NAME}_vs_${NORMAL_BASE_NAME}.cns

# Purpose: Aggregate QC metrics and reports across the run.
module load multiqc
mkdir -p results/multiqc
multiqc results -o results/multiqc
