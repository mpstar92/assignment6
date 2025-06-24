
rule rsem_prepare_reference:
    input:
        "resources/kallisto_transcripts.fasta"
    output:
        directory("results/rsem/index")
    threads: 4
    conda:
        "../envs/rsem.yaml"
    shell:
        """
        rsem-prepare-reference \
            --bowtie2 \
            {input} results/rsem/index/transcripts
        """

rule rsem_calculate_expression:
    input:
        r1 = "results/trimmed/{sample}_R1.fastq.gz",
        r2 = "results/trimmed/{sample}_R2.fastq.gz",
        ref = "results/rsem/index/transcripts"
    output:
        gene = "results/rsem/{sample}.genes.results",
        isoform = "results/rsem/{sample}.isoforms.results"
    log:
        "logs/rsem/{sample}.log"
    threads: 8
    conda:
        "../envs/rsem.yaml"
    shell:
        """
        rsem-calculate-expression \
            --bowtie2 \
            --paired-end {input.r1} {input.r2} \
            {input.ref} results/rsem/{wildcards.sample} \
            --num-threads {threads} &> {log}
        """
