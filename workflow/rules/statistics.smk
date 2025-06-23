rule featurecounts:
    input:
        bam="results/star/{sample}.bam",
        gtf="data/genome/Homo_sapiens.GRCh38.111.gtf"
    output:
        counts="results/counts/{sample}.featureCounts.txt"
    threads: 4
    conda:
        "../envs/subread.yaml"
    shell:
        """
        featureCounts -T {threads} -p -a {input.gtf} -o {output.counts} {input.bam}
        """

#rule htseq_count:
#    input:
#        bam="results/star/{sample}.bam",
#        gtf="data/genome/Homo_sapiens.GRCh38.111.gtf"
#    output:
#        counts="results/counts/{sample}.htseq.txt"
#    threads: 2
#    conda:
#        "../envs/htseq.yaml"
#    shell:
#        """
#        htseq-count -f bam -r pos -s no -t exon -i gene_id {input.bam} {input.gtf} > {output.count}
#        """