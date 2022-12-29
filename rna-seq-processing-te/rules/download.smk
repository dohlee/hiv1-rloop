ASCP_BIN = config['ascp_bin']
ASCP_KEY = config['ascp_key']

rule prefetch_accession:
    output:
        temp('{sample}.sra')
    resources:
        network = 1
    shell:
        f'prefetch --ascp-path "{ASCP_BIN}|{ASCP_KEY}" -v {{wildcards.sample}} --max-size 1000000000 && mv {{wildcards.sample}}/{{wildcards.sample}}.sra . && rm -r {{wildcards.sample}}'

rule parallel_fastq_dump_single:
    input:
        # Required input. Recommend using wildcards for sample names,
        # e.g. {sample,SRR[0-9]+}
        '{sample}.sra'
    output:
        # Required output.
        DATA_DIR / '{sample}.fastq.gz'
    params:
        extra = '--tmpdir .'
    threads: config['threads']['parallel_fastq_dump']
    wrapper:
        'http://dohlee-bio.info:9193/parallel-fastq-dump'

rule parallel_fastq_dump_paired:
    input:
        # Required input. Recommend using wildcards for sample names,
        # e.g. {sample,SRR[0-9]+}
        '{sample}.sra'
    output:
        # Required output.
        DATA_DIR / '{sample}.read1.fastq.gz',
        DATA_DIR / '{sample}.read2.fastq.gz',
    params:
        # Optional parameters. Omit if unused.
        extra = '--tmpdir .'
    threads: config['threads']['parallel_fastq_dump']
    wrapper:
        'http://dohlee-bio.info:9193/parallel-fastq-dump'
