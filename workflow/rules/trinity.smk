rule trinity_assembly:
    input:
        r1 = "results/merged/all_R1.fastq.gz",
        r2 = "results/merged/all_R2.fastq.gz"
    output:
        "results/trinity/Trinity.fasta"
    threads: 8
    conda:
        "../envs/trinity.yaml"
    shell:
        """
        Trinity \
            --seqType fq \
            --left {input.r1} \
            --right {input.r2} \
            --max_memory 20G \
            --CPU {threads} \
            --output results/trinity
        """