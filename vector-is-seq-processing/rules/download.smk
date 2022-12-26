ASCP_BIN = config['ascp_bin']
ASCP_KEY = config['ascp_key']

rule prefetch_accession:
    output:
        temp('{run}.sra')
    resources:
        network = 1
    shell:
        f'prefetch --ascp-path "{ASCP_BIN}|{ASCP_KEY}" -v {{wildcards.run}} && mv {{wildcards.run}}/{{wildcards.run}}.sra . && rm -r {{wildcards.run}}'

rule parallel_fastq_dump_single:
    input:
        # Required input. Recommend using wildcards for sample names,
        # e.g. {sample,SRR[0-9]+}
        '{sample}.sra'
    output:
        # Required output.
        'data/{sample}.fastq.gz'
    params:
        extra = '--tmpdir .'
    threads: 4
    wrapper:
        'http://dohlee-bio.info:9193/parallel-fastq-dump'
