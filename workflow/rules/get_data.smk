rule download_genome:
    output:
        "data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    params:
        url = config["genome"]["fasta"]
    conda:
        "../envs/download.yaml"
    shell:
        """
        wget -O {output} {params.url}
        """
rule download_gtf:
    output:
        "data/genome/Homo_sapiens.GRCh38.111.gtf.gz"
    params:
        url = config["genome"]["gtf"]
    conda:
        "../envs/download.yaml"
    shell:
        """
        wget -O {output} {params.url}
        """

rule unzip_genome:
    input:
        "data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    output:
        "data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    conda:
        "../envs/download.yaml"
    shell:
        """
        gunzip -c {input} > {output}
        """

rule unzip_gtf:
    input:
        "data/genome/Homo_sapiens.GRCh38.111.gtf.gz"
    output:
        "data/genome/Homo_sapiens.GRCh38.111.gtf"
    conda:
        "../envs/download.yaml"
    shell:
        """
        gunzip -c {input} > {output}
        """