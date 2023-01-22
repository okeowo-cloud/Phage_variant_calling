#!/bin/bash

#alias "java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar"

# define samples
SAMPLE="AI-69_S60"
REF_GENOME="T5"
THREADS="8"

# define directories
QUALITY_DIR=data/quality_reports
TRIMMING_DIR=data/trimmed_reads
MAPPED_DIR=mapped
MAPPED_SORTED_DIR=mapped_sorted
MAPPED_STAT_DIR=mapped/mapped_statistics
CALLING_DIR=calling_files

# variables that should be given by
REFERENCE=data/reference/${REF_GENOME}_sequence.fasta
RAW_READ_FOR=data/raw_reads/${SAMPLE}_R1_001.fastq.gz
RAW_READ_REV=data/raw_reads/${SAMPLE}_R2_001.fastq.gz

# variables that we create
TRIMMED_READ_FOR=${TRIMMING_DIR}/trimmed_${SAMPLE}_R1_paired.fq
TRIMMED_READ_REV=${TRIMMING_DIR}/trimmed_${SAMPLE}_R2_paired.fq
TRIMMED_READ_UNPAIRED_FOR=${TRIMMING_DIR}/trimmed_${SAMPLE}_R1_unpaired.fq
TRIMMED_READ_UNPAIRED_REV=${TRIMMING_DIR}/trimmed_${SAMPLE}_R2_inpaired.fq

ALIGNMENT_STAT=${MAPPED_STAT_DIR}/${REF_GENOME}_${SAMPLE}.txt
ALIGNMENT_BAM=${MAPPED_DIR}/${REF_GENOME}_${SAMPLE}.bam
ALIGNMENT_SORTED_BAM=${MAPPED_SORTED_DIR}/${REF_GENOME}_${SAMPLE}_sorted.bam
VARIANTS_VCF=${CALLING_DIR}/variants_ref_${REF_GENOME}_sample_${SAMPLE}_filtered.vcf

# create directories
mkdir ${QUALITY_DIR} ${TRIMMING_DIR}
mkdir ${MAPPED_DIR} ${MAPPED_SORTED_DIR} ${MAPPED_STAT_DIR}
mkdir ${CALLING_DIR}

# quality check
fastqc -q -t ${THREADS} --outdir ${QUALITY_DIR} ${RAW_READ_FOR} ${RAW_READ_REV}

# trimming
java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads ${THREADS} \
    ${RAW_READ_FOR} ${RAW_READ_REV} \
    ${TRIMMED_READ_FOR} ${TRIMMED_READ_UNPAIRED_FOR} \
    ${TRIMMED_READ_REV} ${TRIMMED_READ_UNPAIRED_REV} \
    LEADING:15 TRAILING:15 SLIDINGWINDOW:10:20 MINLEN:20

# index reference
bwa index ${REFERENCE}

# alignment
bwa mem -t ${THREADS} ${REFERENCE} ${TRIMMED_READ_FOR} ${TRIMMED_READ_REV} | \
samtools view -Sb > ${ALIGNMENT_BAM}

# get statistics
samtools flagstat ${ALIGNMENT_BAM} > ${ALIGNMENT_STAT}

# sort alignment
samtools sort -@ ${THREADS} ${ALIGNMENT_BAM} -o ${ALIGNMENT_SORTED_BAM}

# index alignment
samtools index -@ ${THREADS} ${ALIGNMENT_SORTED_BAM} ${ALIGNMENT_SORTED_BAM}.bai

# variant calling
bcftools mpileup -Ou -f ${REFERENCE} ${ALIGNMENT_SORTED_BAM} | \
bcftools call -Ou -mv --ploidy 1 | \
bcftools filter -s LowQual > ${VARIANTS_VCF}
