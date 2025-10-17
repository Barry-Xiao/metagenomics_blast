from snakemake.io import glob_wildcards

# load configuration file
configfile: "config.yaml"

# detect gzipped fastq files in the input directory
SAMPLES, = glob_wildcards(config["input_folder"] + "/{sample}.fastq.gz")

# helper function for conditional temp()
def temp_decider(file_path,
                 keep_intermediate=config.get("keep_intermediate", False)
                ):
    if keep_intermediate:
        return file_path
    else:
        return temp(file_path)

# main pipeline

## rule all: target rule, annotated tsv with taxonomy
rule all:
    input:
        "results/merged_annotated_hits.tsv"

## rule subsample_fastq: subsample the input fastq files
rule subsample_fastq:
    input:
        config["input_folder"] + "/{sample}.fastq.gz"
    output:
        temp_decider("results/subsample/{sample}.fastq.gz")
    params:
        num_reads=config["subsample_reads"]
    log:
        "logs/subsample/{sample}.log"
    conda:
        "envs/subsample_blast.yaml"
    shell:
        """
        seqtk sample -s 20251016 {input} {params.num_reads} | gzip > {output} 2> {log}
        """

## rule fastq_to_fasta: convert fastq to fasta
rule fastq_to_fasta:
    input:
        "results/subsample/{sample}.fastq.gz"
    output:
        temp_decider("results/fasta/{sample}.fasta")
    log:
        "logs/fasta/{sample}.log"
    conda:
        "envs/subsample_blast.yaml"
    shell:
        """
        seqtk seq -a {input} > {output} 2> {log}
        """

## rule blastn: aligns sequences using blastn against the database
rule blastn:
    input:
        fasta="results/fasta/{sample}.fasta"
    output:
        temp_decider("results/blastn/{sample}.blastn.tsv")
    log:
        "logs/blastn/{sample}.log"
    conda:
        "envs/subsample_blast.yaml"
    shell:
        """
        blastn -query {input.fasta} \\
               -db {config[remote_db]} \\
               -out {output} \\
               -outfmt 6 \\
               -remote \\
               2> {log}
        """

## rule find best hit: parse the BLAST output to find the best hit
rule find_best_hit:
    input:
        blast_result="results/blastn/{sample}.blastn.tsv"
    output:
        temp_decider("results/best_hits/{sample}.best_hit.tsv")
    log:
        "logs/find_best_hit/{sample}.log"
    conda:
        "envs/subsample_blast.yaml"
    script:
        "scripts/find_best_hit.py"


## rule get taxonomy: get taxonomy for the matched sseqid
rule get_taxonomy:
    input:
        "results/best_hits/{sample}.best_hit.tsv"
    output:
        temp_decider("results/taxonomy/{sample}.taxonomy.tsv")
    log:
        "logs/get_taxonomy/{sample}.log"
    conda:
        "envs/subsample_blast.yaml"
    retries: 3  # Automatically retry the job up to 3 times upon failure.
    run:
        import os
        if os.path.getsize(input[0]) < 100:
            shell("touch {output}")
        else:
            shell(
                """
                timeout 1800 cut -f 2 {input} | \\
                tail -n +2 | \\
                sort -u | \\
                epost -db nuccore | \\
                efetch -format docsum | \\
                xtract -pattern DocumentSummary -element AccessionVersion,Organism,TaxId \\
                > {output} 2> {log}
                """
            )

## rule annotate with taxonomy: merge the best hit results with fetched taxonomy info
rule annotate_with_taxonomy:
    input:
        hits="results/best_hits/{sample}.best_hit.tsv",
        tax="results/taxonomy/{sample}.taxonomy.tsv"
    output:
        temp_decider("results/annotated_hits/{sample}.annotated.tsv")
    log:
        "logs/annotate/{sample}.log"
    conda:
        "envs/subsample_blast.yaml"
    script:
        "scripts/merge_taxonomy.py"

## rule merge outputs: merge all outputs for all samples together
rule merge_outputs:
    input:
        expand("results/annotated_hits/{sample}.annotated.tsv", sample=SAMPLES)
    output:
        "results/merged_annotated_hits.tsv"
    log:
        "logs/merge_outputs.log"
    conda:
        "envs/subsample_blast.yaml"
    script:
        "scripts/merge_outputs.py"


