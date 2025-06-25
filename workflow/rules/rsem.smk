rule rsem_prepare_reference_trinity:
    input:
        "results/trinity/Trinity.fasta"
    output:
        "results/rsem_trinity/index/Trinity.grp"
    log:
        "logs/rsem_trinity/index.log"
    threads: 4
    conda:
        "../envs/rsem.yaml"
    shell:
        """
        rsem-prepare-reference \
            --bowtie2 \
            {input} results/rsem_trinity/index/Trinity &> {log}
        """

rule rsem_quant_trinity:
    input:
        r1 = "results/merged/all_R1.fastq.gz",
        r2 = "results/merged/all_R2.fastq.gz",
        ref = "results/rsem_trinity/index/Trinity.grp"
    output:
        gene = "results/rsem_trinity/quant.genes.results",
        isoform = "results/rsem_trinity/quant.isoforms.results"
    log:
        "logs/rsem_trinity/quant.log"
    threads: 8
    conda:
        "../envs/rsem.yaml"
    shell:
        """
        rsem-calculate-expression \
            --bowtie2 \
            --paired-end {input.r1} {input.r2} \
            {input.ref} results/rsem_trinity/quant \
            --num-threads {threads} &> {log}
        """