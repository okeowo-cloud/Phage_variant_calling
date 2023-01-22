import os

rule all:
    input:
        expand("calling_files/{VCF_FILE}", VCF_FILE="variants_ref_T5_sample_AI-72_S63_filtered.vcf"),
        expand("data/quality_reports/{QUALITY_FILES}", QUALITY_FILES="AI-72_S63_R1_001_fastqc.html" ),
        expand("data/trimmed_reads/trimmed_AI-72_S63_{TRIMMED_FILES}",
            TRIMMED_FILES=["R1_paired.fq", "R1_unpaired.fq", "R2_paired.fq", "R2_unpaired.fq"]),
        expand("data/reference/{REFERENCE_FILE}", REFERENCE_FILE="T5_sequence.fasta.ann"),
        os.path.join("mapped/mapped_statistics", "T5_AI-72_S63.txt"),
        os.path.join("mapped/mapped_statistics", "T5_AI-72_S63_sorted.bam.bai")

rule create_dir:
    output:
        directory(expand("{dirs}",
            dirs=["mapped/mapped_statistics", "data/quality_reports", "data/trimmed_reads", "mapped_sorted", "calling_files"]))
    shell:
        "mkdir {output}"


rule check_read_quality:
    params:
        output_dir="data/quality_reports"
    input:
        os.path.join("data/raw_reads", "AI-72_S63_R1_001.fastq.gz"),
        os.path.join("data/raw_reads", "AI-72_S63_R2_001.fastq.gz")
    output:
        os.path.join("data/quality_reports","{QUALITY_FILES}")
    threads: 8
    conda: os.path.join("envs", "fastqc_env.yml")
    shell:
        "fastqc -q -t {threads} --outdir {params.output_dir} {input}"


rule trim_reads:
    params:
        leading = 15,
        trailing = 15,
        length = 10,
        breadth = 20,
        minlen = 20
    input:
        os.path.join("data/raw_reads","AI-72_S63_R1_001.fastq.gz"),
        os.path.join("data/raw_reads","AI-72_S63_R2_001.fastq.gz")
    output:
        os.path.join("data/trimmed_reads" ,"trimmed_AI-72_S63_{TRIMMED_FILES}")
    conda: os.path.join("envs", "trimmomatic_env.yml")
    threads: 8
    shell:
        "java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads {threads} {input} {output} \
        LEADING:{params.leading} \
        TRAILING:{params.trailing} \
        SLIDINGWINDOW:{params.length}:{params.breadth} \
        MINLEN:{params.minlen}"


rule index_ref:
    input:
        os.path.join("data/reference", "T5_sequence.fasta")
    output:
        os.path.join("data/reference", "T5_sequence.fasta.ann")
    threads: 8
    conda: os.path.join("envs", "bwa_env.yml")
    shell:
        "bwa index {input}"


rule alignment:
    input:
        os.path.join("data/reference", "T5_sequence.fasta"),
        os.path.join("data/raw_reads","AI-72_S63_R1_001.fastq.gz"),
        os.path.join("data/raw_reads","AI-72_S63_R2_001.fastq.gz")
    output:
        os.path.join("mapped", "T5_AI-72_S63.bam")
    threads: 8
    conda: os.path.join("envs", "bwa_env.yml")
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb > {output}"


rule get_statistics:
    input:
        os.path.join("mapped","T5_AI-72_S63.bam")
    output:
        os.path.join("mapped/mapped_statistics", "T5_AI-72_S63.txt")
    conda: os.path.join("envs", "samtools_env.yml")
    shell:
        "samtools flagstat {input} > {output}"

rule sort_alignment:
    input:
        os.path.join("mapped","T5_AI-72_S63.bam")
    output:
        os.path.join("mapped_sorted", "T5_AI-72_S63_sorted.bam")
    conda: os.path.join("envs", "samtools_env.yml")
    threads: 8
    shell:
        "samtools sort -@ {threads} {input} -o {output}"


rule index_alignment:
    input:
        os.path.join("mapped_sorted", "T5_AI-72_S63_sorted.bam")
    output:
        os.path.join("mapped/mapped_statistics", "T5_AI-72_S63_sorted.bam.bai")
    threads: 8
    conda: os.path.join("envs", "samtools_env.yml")
    shell:
        "samtools index -@ {threads} {input} {output}"

rule variant_calling:
    input:
        os.path.join("data/reference", "T5_sequence.fasta"),
        os.path.join("mapped_sorted","T5_AI-72_S63_sorted.bam")
    output:
        os.path.join("calling_files", "{VCF_FILE}")
    conda: os.path.join("envs", "bcftools_env.yml")
    shell:
        "bcftools mpileup -Ou -f {input} | bcftools call -Ou -mv --ploidy 1 | bcftools filter -s LowQual > {output}"
