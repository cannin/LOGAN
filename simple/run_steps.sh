#!/usr/bin/env bash
set -euo pipefail

# Run from the simple/ directory so relative paths resolve.

mkdir -p results/fastp
fastp -w 8 \
    -i /data/fastq/TUMOR_1_R1.fastq.gz -I /data/fastq/TUMOR_1_R2.fastq.gz \
    -o results/fastp/TUMOR_1.R1.trimmed.fastq.gz -O results/fastp/TUMOR_1.R2.trimmed.fastq.gz \
    --json results/fastp/TUMOR_1.fastp.json \
    --html results/fastp/TUMOR_1.fastp.html
fastp -w 8 \
    -i /data/fastq/NORMAL_1_R1.fastq.gz -I /data/fastq/NORMAL_1_R2.fastq.gz \
    -o results/fastp/NORMAL_1.R1.trimmed.fastq.gz -O results/fastp/NORMAL_1.R2.trimmed.fastq.gz \
    --json results/fastp/NORMAL_1.fastp.json \
    --html results/fastp/NORMAL_1.fastp.html

mkdir -p results/qc/fastqc
fastqc -t 8 -o results/qc/fastqc \
    results/fastp/TUMOR_1.R1.trimmed.fastq.gz results/fastp/TUMOR_1.R2.trimmed.fastq.gz
fastqc -t 8 -o results/qc/fastqc \
    results/fastp/NORMAL_1.R1.trimmed.fastq.gz results/fastp/NORMAL_1.R2.trimmed.fastq.gz

mkdir -p results/qc/fastq_screen
fastq_screen --conf ../conf/fastq_screen.conf \
    --outdir results/qc/fastq_screen \
    --threads 8 \
    --subset 1000000 \
    --aligner bowtie2 \
    --force \
    results/fastp/TUMOR_1.R1.trimmed.fastq.gz results/fastp/TUMOR_1.R2.trimmed.fastq.gz
fastq_screen --conf ../conf/fastq_screen.conf \
    --outdir results/qc/fastq_screen \
    --threads 8 \
    --subset 1000000 \
    --aligner bowtie2 \
    --force \
    results/fastp/NORMAL_1.R1.trimmed.fastq.gz results/fastp/NORMAL_1.R2.trimmed.fastq.gz

mkdir -p results/qc/kraken
kraken2 --db /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/kraken/20180907_standard_kraken2 \
    --threads 8 --report results/qc/kraken/TUMOR_1.trimmed.kraken_bacteria.taxa.txt \
    --output - \
    --gzip-compressed \
    --paired results/fastp/TUMOR_1.R1.trimmed.fastq.gz results/fastp/TUMOR_1.R2.trimmed.fastq.gz
cut -f2,3 results/qc/kraken/TUMOR_1.trimmed.kraken_bacteria.taxa.txt | \
    ktImportTaxonomy - -o results/qc/kraken/TUMOR_1.trimmed.kraken_bacteria.krona.html
kraken2 --db /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/kraken/20180907_standard_kraken2 \
    --threads 8 --report results/qc/kraken/NORMAL_1.trimmed.kraken_bacteria.taxa.txt \
    --output - \
    --gzip-compressed \
    --paired results/fastp/NORMAL_1.R1.trimmed.fastq.gz results/fastp/NORMAL_1.R2.trimmed.fastq.gz
cut -f2,3 results/qc/kraken/NORMAL_1.trimmed.kraken_bacteria.taxa.txt | \
    ktImportTaxonomy - -o results/qc/kraken/NORMAL_1.trimmed.kraken_bacteria.krona.html

mkdir -p results/align
mkdir -p /tmp/logan_simple/bwa_TUMOR_1
if lscpu | grep -q avx2; then
    bwa-mem2.avx2 mem -M \
        -R '@RG\tID:TUMOR_1\tSM:TUMOR_1\tPL:illumina\tLB:TUMOR_1\tPU:TUMOR_1\tDS:wgs' \
        -t 8 \
        /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
        results/fastp/TUMOR_1.R1.trimmed.fastq.gz results/fastp/TUMOR_1.R2.trimmed.fastq.gz | \
    samblaster -M | \
    samtools sort -T /tmp/logan_simple/bwa_TUMOR_1/ \
        -@ 8 -m 2G - \
        --write-index -o results/align/TUMOR_1.bam##idx##results/align/TUMOR_1.bam.bai
elif lscpu | grep -q sse4_1; then
    bwa-mem2.sse41 mem -M \
        -R '@RG\tID:TUMOR_1\tSM:TUMOR_1\tPL:illumina\tLB:TUMOR_1\tPU:TUMOR_1\tDS:wgs' \
        -t 8 \
        /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
        results/fastp/TUMOR_1.R1.trimmed.fastq.gz results/fastp/TUMOR_1.R2.trimmed.fastq.gz | \
    samblaster -M | \
    samtools sort -T /tmp/logan_simple/bwa_TUMOR_1/ \
        -@ 8 -m 2G - \
        --write-index -o results/align/TUMOR_1.bam##idx##results/align/TUMOR_1.bam.bai
else
    bwa-mem2 mem -M \
        -R '@RG\tID:TUMOR_1\tSM:TUMOR_1\tPL:illumina\tLB:TUMOR_1\tPU:TUMOR_1\tDS:wgs' \
        -t 8 \
        /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
        results/fastp/TUMOR_1.R1.trimmed.fastq.gz results/fastp/TUMOR_1.R2.trimmed.fastq.gz | \
    samblaster -M | \
    samtools sort -T /tmp/logan_simple/bwa_TUMOR_1/ \
        -@ 8 -m 2G - \
        --write-index -o results/align/TUMOR_1.bam##idx##results/align/TUMOR_1.bam.bai
fi

mkdir -p results/align
mkdir -p /tmp/logan_simple/bwa_NORMAL_1
if lscpu | grep -q avx2; then
    bwa-mem2.avx2 mem -M \
        -R '@RG\tID:NORMAL_1\tSM:NORMAL_1\tPL:illumina\tLB:NORMAL_1\tPU:NORMAL_1\tDS:wgs' \
        -t 8 \
        /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
        results/fastp/NORMAL_1.R1.trimmed.fastq.gz results/fastp/NORMAL_1.R2.trimmed.fastq.gz | \
    samblaster -M | \
    samtools sort -T /tmp/logan_simple/bwa_NORMAL_1/ \
        -@ 8 -m 2G - \
        --write-index -o results/align/NORMAL_1.bam##idx##results/align/NORMAL_1.bam.bai
elif lscpu | grep -q sse4_1; then
    bwa-mem2.sse41 mem -M \
        -R '@RG\tID:NORMAL_1\tSM:NORMAL_1\tPL:illumina\tLB:NORMAL_1\tPU:NORMAL_1\tDS:wgs' \
        -t 8 \
        /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
        results/fastp/NORMAL_1.R1.trimmed.fastq.gz results/fastp/NORMAL_1.R2.trimmed.fastq.gz | \
    samblaster -M | \
    samtools sort -T /tmp/logan_simple/bwa_NORMAL_1/ \
        -@ 8 -m 2G - \
        --write-index -o results/align/NORMAL_1.bam##idx##results/align/NORMAL_1.bam.bai
else
    bwa-mem2 mem -M \
        -R '@RG\tID:NORMAL_1\tSM:NORMAL_1\tPL:illumina\tLB:NORMAL_1\tPU:NORMAL_1\tDS:wgs' \
        -t 8 \
        /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
        results/fastp/NORMAL_1.R1.trimmed.fastq.gz results/fastp/NORMAL_1.R2.trimmed.fastq.gz | \
    samblaster -M | \
    samtools sort -T /tmp/logan_simple/bwa_NORMAL_1/ \
        -@ 8 -m 2G - \
        --write-index -o results/align/NORMAL_1.bam##idx##results/align/NORMAL_1.bam.bai
fi

mkdir -p results/bqsr
gatk --java-options '-Xmx10g' BaseRecalibrator \
    --input results/align/TUMOR_1.bam \
    --reference /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    --known-sites /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/GATK_resource_bundle/dbsnp_138.hg38.vcf.gz \
    --known-sites /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-sites /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/GATK_resource_bundle/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz \
    --intervals /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list \
    --output results/bqsr/TUMOR_1.recal_data.grp
gatk --java-options '-Xmx10g' BaseRecalibrator \
    --input results/align/NORMAL_1.bam \
    --reference /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    --known-sites /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/GATK_resource_bundle/dbsnp_138.hg38.vcf.gz \
    --known-sites /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-sites /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/GATK_resource_bundle/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz \
    --intervals /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list \
    --output results/bqsr/NORMAL_1.recal_data.grp

gatk --java-options '-Xmx20g' ApplyBQSR \
    --reference /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    --input results/align/TUMOR_1.bam \
    --bqsr-recal-file results/bqsr/TUMOR_1.recal_data.grp \
    --output results/align/TUMOR_1.bqsr.bam \
    --use-jdk-inflater \
    --use-jdk-deflater
samtools index -@ 8 results/align/TUMOR_1.bqsr.bam results/align/TUMOR_1.bqsr.bam.bai
gatk --java-options '-Xmx20g' ApplyBQSR \
    --reference /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    --input results/align/NORMAL_1.bam \
    --bqsr-recal-file results/bqsr/NORMAL_1.recal_data.grp \
    --output results/align/NORMAL_1.bqsr.bam \
    --use-jdk-inflater \
    --use-jdk-deflater
samtools index -@ 8 results/align/NORMAL_1.bqsr.bam results/align/NORMAL_1.bqsr.bam.bai

mkdir -p results/qc/align
samtools flagstat results/align/TUMOR_1.bqsr.bam > results/qc/align/TUMOR_1.samtools_flagstat.txt
samtools flagstat results/align/NORMAL_1.bqsr.bam > results/qc/align/NORMAL_1.samtools_flagstat.txt

mkdir -p results/qc/align
mosdepth -n --fast-mode --by 500 results/qc/align/TUMOR_1 results/align/TUMOR_1.bqsr.bam -t 8
mv results/qc/align/TUMOR_1.mosdepth.region.dist.txt results/qc/align/TUMOR_1.mosdepth.region.dist.txt
mv results/qc/align/TUMOR_1.mosdepth.summary.txt results/qc/align/TUMOR_1.mosdepth.summary.txt
mv results/qc/align/TUMOR_1.regions.bed.gz results/qc/align/TUMOR_1.regions.bed.gz
mv results/qc/align/TUMOR_1.regions.bed.gz.csi results/qc/align/TUMOR_1.regions.bed.gz.csi
mosdepth -n --fast-mode --by 500 results/qc/align/NORMAL_1 results/align/NORMAL_1.bqsr.bam -t 8
mv results/qc/align/NORMAL_1.mosdepth.region.dist.txt results/qc/align/NORMAL_1.mosdepth.region.dist.txt
mv results/qc/align/NORMAL_1.mosdepth.summary.txt results/qc/align/NORMAL_1.mosdepth.summary.txt
mv results/qc/align/NORMAL_1.regions.bed.gz results/qc/align/NORMAL_1.regions.bed.gz
mv results/qc/align/NORMAL_1.regions.bed.gz.csi results/qc/align/NORMAL_1.regions.bed.gz.csi

mkdir -p results/qc/align
unset DISPLAY
qualimap bamqc -bam results/align/TUMOR_1.bqsr.bam \
    --java-mem-size=96G \
    -c -ip \
    -outdir results/qc/align/TUMOR_1 \
    -outformat HTML \
    -nt 8 \
    -nw 500 \
    -p NON-STRAND-SPECIFIC
mv results/qc/align/TUMOR_1/genome_results.txt results/qc/align/TUMOR_1_genome_results.txt
mv results/qc/align/TUMOR_1/qualimapReport.html results/qc/align/TUMOR_1_qualimapReport.html
qualimap bamqc -bam results/align/NORMAL_1.bqsr.bam \
    --java-mem-size=96G \
    -c -ip \
    -outdir results/qc/align/NORMAL_1 \
    -outformat HTML \
    -nt 8 \
    -nw 500 \
    -p NON-STRAND-SPECIFIC
mv results/qc/align/NORMAL_1/genome_results.txt results/qc/align/NORMAL_1_genome_results.txt
mv results/qc/align/NORMAL_1/qualimapReport.html results/qc/align/NORMAL_1_qualimapReport.html

mkdir -p results/snv/TUMOR_1_vs_NORMAL_1/mutect2
gatk Mutect2 \
    --reference /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    --intervals /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list \
    --input results/align/TUMOR_1.bqsr.bam \
    --input results/align/NORMAL_1.bqsr.bam \
    --normal-sample NORMAL_1 \
    --tumor-sample TUMOR_1 \
    --germline-resource /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/GATK_GRCh38.d1.vd1/somatic-hg38-af-only-gnomad.hg38.vcf.gz \
    --panel-of-normals /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/PON/MuTect2.PON.5210.vcf.gz \
    --output results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.mutect2.vcf.gz \
    --f1r2-tar-gz results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.f1r2.tar.gz \
    --independent-mates

mkdir -p results/snv/TUMOR_1_vs_NORMAL_1/mutect2
gatk --java-options -Xmx48g GetPileupSummaries \
    -I results/align/TUMOR_1.bqsr.bam \
    -V /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/GATK_GRCh38.d1.vd1/somatic-hg38-af-only-gnomad.hg38.vcf.gz \
    --intervals /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list \
    -O results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.tpileup.table
gatk --java-options -Xmx48g GetPileupSummaries \
    -I results/align/NORMAL_1.bqsr.bam \
    -V /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/GATK_GRCh38.d1.vd1/somatic-hg38-af-only-gnomad.hg38.vcf.gz \
    --intervals /data/CCBR_Pipeliner/Pipelines/LOGAN/resources/hg38/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list \
    -O results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.npileup.table

mkdir -p results/snv/TUMOR_1_vs_NORMAL_1/mutect2
gatk GatherPileupSummaries \
    --sequence-dictionary /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem/GRCh38.d1.vd1.dict \
    -I results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.tpileup.table \
    -O results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.allpileups.table
gatk GatherPileupSummaries \
    --sequence-dictionary /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem/GRCh38.d1.vd1.dict \
    -I results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.npileup.table \
    -O results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.normal.allpileups.table
gatk CalculateContamination \
    -I results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.allpileups.table \
    --matched-normal results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.normal.allpileups.table \
    -O results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.contamination.table
gatk CalculateContamination \
    -I results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.normal.allpileups.table \
    -O results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.normal.contamination.table

mkdir -p results/snv/TUMOR_1_vs_NORMAL_1/mutect2
gatk LearnReadOrientationModel \
    --output results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.read-orientation-model.tar.gz \
    --input results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.f1r2.tar.gz

mkdir -p results/snv/TUMOR_1_vs_NORMAL_1/mutect2
gatk MergeMutectStats \
    --stats results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.mutect2.vcf.gz.stats \
    -O results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.final.stats

mkdir -p results/snv/TUMOR_1_vs_NORMAL_1/mutect2
gatk SortVcf \
    -I results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.mutect2.vcf.gz \
    -O results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.mutect2.concat.vcf.gz \
    --CREATE_INDEX
gatk FilterMutectCalls \
    -R /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    -V results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.mutect2.concat.vcf.gz \
    --ob-priors results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.read-orientation-model.tar.gz \
    --contamination-table results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.contamination.table \
    --stats results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.final.stats \
    -O results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.mutect2.marked.vcf.gz
gatk SelectVariants \
    -R /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    --variant results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.mutect2.marked.vcf.gz \
    --exclude-filtered \
    --output results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.mutect2.final.vcf.gz
bcftools sort results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.mutect2.final.vcf.gz | \
    bcftools norm --threads 8 --check-ref s \
    -f /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa -O v | \
    awk '{gsub(/\y[W|K|Y|R|S|M|B|D|H|V]\y/,"N",$4); OFS = "\t"; print}' | \
    sed '/^$/d' | \
    bcftools view - -Oz -o results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.mutect2.norm.vcf.gz
bcftools index -t results/snv/TUMOR_1_vs_NORMAL_1/mutect2/TUMOR_1_vs_NORMAL_1.mutect2.norm.vcf.gz

mkdir -p results/sv/TUMOR_1_vs_NORMAL_1/manta/run
configManta.py \
    --tumorBam results/align/TUMOR_1.bqsr.bam \
    --normalBam results/align/NORMAL_1.bqsr.bam \
    --referenceFasta /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    --runDir results/sv/TUMOR_1_vs_NORMAL_1/manta/run
results/sv/TUMOR_1_vs_NORMAL_1/manta/run/runWorkflow.py -m local -j 8
cp results/sv/TUMOR_1_vs_NORMAL_1/manta/run/results/variants/somaticSV.vcf.gz \
    results/sv/TUMOR_1_vs_NORMAL_1/manta/TUMOR_1_vs_NORMAL_1.somaticSV.vcf.gz
cp results/sv/TUMOR_1_vs_NORMAL_1/manta/run/results/variants/somaticSV.vcf.gz.tbi \
    results/sv/TUMOR_1_vs_NORMAL_1/manta/TUMOR_1_vs_NORMAL_1.somaticSV.vcf.gz.tbi

mkdir -p results/cnv/TUMOR_1_vs_NORMAL_1/cnvkit
cnvkit.py batch \
    results/align/TUMOR_1.bqsr.bam -n results/align/NORMAL_1.bqsr.bam \
    -f /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/bwamem2/GRCh38.d1.vd1.fa \
    -m wgs \
    --access /data/CCBR_Pipeliner/Pipelines/XAVIER/resources/hg38/cnvkit/access-10kb.hg38.bed \
    --output-reference results/cnv/TUMOR_1_vs_NORMAL_1/cnvkit/TUMOR_1_vs_NORMAL_1.cnn \
    --output-dir results/cnv/TUMOR_1_vs_NORMAL_1/cnvkit
cp results/cnv/TUMOR_1_vs_NORMAL_1/cnvkit/TUMOR_1.bqsr.cns \
    results/cnv/TUMOR_1_vs_NORMAL_1/cnvkit/TUMOR_1_vs_NORMAL_1.cns

mkdir -p results/multiqc
multiqc results -o results/multiqc
