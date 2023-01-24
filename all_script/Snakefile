import os

"""
Snakemake version 6.10.0 <specify your version here>
Python version 3.10.9 <specify your version here>

To run this code, write in terminal:
snakemake --use-conda -c <specify number of threads> -n 
"""

# define samples
SAMPLES = ['AI-69_S60', 'AI-70_S61', 'AI-71_S62', 'AI-72_S63', 'AI-73_S64']
REFERENCES = ['T5']

# define directories
reference_dir = os.path.join('data', 'reference')
raw_reads_dir = os.path.join('data', 'raw_reads')
mapped_dir = os.path.join('mapped')
sorted_dir = os.path.join('mapped_sorted')
calling_dir = os.path.join('calling')
result_dir = os.path.join('results')

# create directories
os.makedirs(reference_dir, exist_ok=True)
os.makedirs(raw_reads_dir, exist_ok=True)
os.makedirs(mapped_dir, exist_ok=True)
os.makedirs(sorted_dir, exist_ok=True)
os.makedirs(calling_dir, exist_ok=True)
os.makedirs(result_dir, exist_ok=True)


rule all:
    input:
        expand(os.path.join(result_dir, 'sorted_{reference}_sample_{sample}_flagstat.log'), reference=REFERENCES, sample=SAMPLES),
        expand(os.path.join(calling_dir, 'variants_filtered_{reference}_sample_{sample}.vcf'), reference=REFERENCES, sample=SAMPLES)

rule index_reference:
    input:
        os.path.join(reference_dir, '{reference}_sequence.fasta')
    output:
        os.path.join(reference_dir, '{reference}_sequence.fasta.amb'),
        os.path.join(reference_dir, '{reference}_sequence.fasta.ann'),
        os.path.join(reference_dir, '{reference}_sequence.fasta.bwt'),
        os.path.join(reference_dir, '{reference}_sequence.fasta.pac'),
        os.path.join(reference_dir, '{reference}_sequence.fasta.sa')
    conda:
        os.path.join('envs', 'bwa_env.yml')
    shell:
        "bwa index {input}"

rule bwa_map:
    input:
        index_amb = os.path.join(reference_dir, '{reference}_sequence.fasta.amb'),
        index_ann = os.path.join(reference_dir, '{reference}_sequence.fasta.ann'),
        index_bwt = os.path.join(reference_dir, '{reference}_sequence.fasta.bwt'),
        index_pac = os.path.join(reference_dir, '{reference}_sequence.fasta.pac'),
        index_sa = os.path.join(reference_dir, '{reference}_sequence.fasta.sa'),
        fasta = os.path.join(reference_dir, '{reference}_sequence.fasta'),
        fastq1 = os.path.join(raw_reads_dir, '{sample}_R1_001.fastq.gz'),
        fastq2 = os.path.join(raw_reads_dir, '{sample}_R2_001.fastq.gz')
    output:
        bam = os.path.join(mapped_dir, '{reference}_sample_{sample}.bam')
    conda:
        os.path.join('envs', 'bwa_env.yml')
    threads: 8
    shell:
        "bwa mem -t {threads} {input.fasta} {input.fastq1} {input.fastq2} | samtools view -Sb > {output.bam}"

rule get_statistics:
    input:
        os.path.join(sorted_dir, 'sorted_{reference}_sample_{sample}.bam')
    output:
        os.path.join(result_dir, 'sorted_{reference}_sample_{sample}_flagstat.log')
    threads: 1
    conda:
        os.path.join('envs', 'samtools_env.yml')
    shell:
        "samtools flagstat {input} > {output}"

rule sam_sort:
    input:
        os.path.join(mapped_dir, '{reference}_sample_{sample}.bam')
    output:
        os.path.join(sorted_dir, 'sorted_{reference}_sample_{sample}.bam')
    threads: 7
    conda:
        os.path.join('envs', 'samtools_env.yml')
    shell:
        "samtools sort -@ {threads} {input} -o {output}"

rule sam_index:
    input:
        os.path.join(sorted_dir, 'sorted_{reference}_sample_{sample}.bam')
    output:
        os.path.join(sorted_dir, 'sorted_{reference}_sample_{sample}.bam.bai')
    threads: 8
    conda:
        os.path.join('envs', 'samtools_env.yml')
    shell:
        "samtools index -@ {threads} {input} {output}"

rule variant_calling:
    input:
        fasta=os.path.join(reference_dir, '{reference}_sequence.fasta'),
        bam=os.path.join(sorted_dir, 'sorted_{reference}_sample_{sample}.bam'),
        index_bai=os.path.join(sorted_dir, 'sorted_{reference}_sample_{sample}.bam.bai')
    output:
        os.path.join(calling_dir, 'variants_filtered_{reference}_sample_{sample}.vcf')
    conda:
        os.path.join('envs', 'bcftools_env.yml')
    shell:
        """
        bcftools mpileup -Ou -f {input.fasta} {input.bam} | \
        bcftools call -Ou -mv --ploidy 1 | \
        bcftools filter -s LowQual > {output}
        """
