rule star_index:
    input:
        fasta="data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        gtf="data/genome/Homo_sapiens.GRCh38.111.gtf"
    output:
        directory("data/genome/star_index")
    threads: 8
    conda:
        "../envs/star.yaml"
    shell:
        """
        STAR --runThreadN {threads} --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf}
        """

rule star_align:
    input:
        index="data/genome/star_index",
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz"
    output:
        bam=temp("results/star/{sample}.bam")
    threads: 8
    conda:
        "../envs/star.yaml"
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.index} \
             --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat \
             --outSAMtype BAM Unsorted --outFileNamePrefix results/star/{wildcards.sample}_ \
        && mv results/star/{wildcards.sample}_Aligned.out.bam {output.bam}
        """