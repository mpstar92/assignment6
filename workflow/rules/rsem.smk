rule rsem_prepare_reference_trinity:
    input:
        assembly = "results/trinity.Trinity.fasta",
        gene_trans_map = "results/trinity.Trinity.fasta.gene_trans_map"
    output:
        grp = "results/rsem_trinity/index/Trinity.grp",
        bt2_idx = expand("results/rsem_trinity/index/Trinity.{ext}.bt2", 
                        ext=["1", "2", "3", "4", "rev.1", "rev.2"])
    threads: 8
    conda:
        "../envs/rsem.yaml"
    shell:
        """
        rsem-prepare-reference --transcript-to-gene-map {input.gene_trans_map} --bowtie2 {input.assembly} results/rsem_trinity/index/Trinity
        """

rule rsem_quant_trinity:
    input:
        r1 = "results/merged/all_R1.fastq.gz",
        r2 = "results/merged/all_R2.fastq.gz",
        grp = "results/rsem_trinity/index/Trinity.grp",
        bt2_idx = expand("results/rsem_trinity/index/Trinity.{ext}.bt2", 
                        ext=["1", "2", "3", "4", "rev.1", "rev.2"])
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
            results/rsem_trinity/index/Trinity results/rsem_trinity/quant \
            --num-threads {threads} &> {log}
        """