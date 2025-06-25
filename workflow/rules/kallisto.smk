rule kallisto_index_trinity:
    input:
        "results/trinity/Trinity.fasta"
    output:
        "results/kallisto_trinity/index/Trinity.idx"
    threads: 1
    log:
        "logs/kallisto_trinity/index.log"
    conda:
        "../envs/kallisto.yaml"
    shell:
        """
        kallisto index -i {output} {input} &> {log}
        """

rule kallisto_quant_trinity:
    input:
        index = "results/kallisto_trinity/index/Trinity.idx",
        r1 = "results/merged/all_R1.fastq.gz",
        r2 = "results/merged/all_R2.fastq.gz"
    output:
        "results/kallisto_trinity/quant/abundance.tsv"
    threads: 20
    log:
        "logs/kallisto_trinity/quant.log"
    conda:
        "../envs/kallisto.yaml"
    shell:
        """
        kallisto quant \
            -i {input.index} \
            -o {output} \
            -b 100 -t {threads} \
            {input.r1} {input.r2} &> {log}
        """