
rule kallisto_index_provided:
    input:
        "resources/kallisto_transcripts.fasta"
    output:
        "results/kallisto/index/transcripts.idx"
    threads: 1
    conda:
        "../envs/kallisto.yaml"
    shell:
        """
        kallisto index -i {output} {input}
        """

rule kallisto_quant:
    input:
        index = "results/kallisto/index/transcripts.idx",
        r1 = "results/trimmed/{sample}_R1.fastq.gz",
        r2 = "results/trimmed/{sample}_R2.fastq.gz"
    output:
        abundance = "results/kallisto/{sample}/abundance.tsv"
    log:
        "logs/kallisto/{sample}.log"
    threads: 8
    conda:
        "../envs/kallisto.yaml"
    shell:
        """
        kallisto quant \
            -i {input.index} \
            -o results/kallisto/{wildcards.sample}/ \
            -b 100 \
            -t {threads} \
            {input.r1} {input.r2} &> {log}
        """