rule merge_reads:
    input:
        r1 = expand("results/trimmed/{sample}_R1.fastq.gz", sample=SAMPLES),
        r2 = expand("results/trimmed/{sample}_R2.fastq.gz", sample=SAMPLES)
    output:
        r1 = "results/merged/all_R1.fastq.gz",
        r2 = "results/merged/all_R2.fastq.gz"
    conda:
        "../envs/base.yaml"
    shell:
        """
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
        """