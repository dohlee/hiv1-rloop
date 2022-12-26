
rule macs2_filterdup:
    input:
        result_dir / '02_bwa' / '{sample}.sorted.bam'
    output:
        result_dir / '03_macs2' / '{sample}.sorted.filterdup.bed'
    params:
        # Extra options.
        extra = '',
        # Mappable genome size. (Available preset: hs, mm, ce, dm)
        # Default: hs
        gsize = 'hs',
        # Tag size. This will override the auto detected tag size.
        # Default: Not set
        tsize = False,
        # Pvalue cutoff for binomial distribution test.
        # Default: 1e-5
        pvalue = 1e-5,
        # MACS2 filterdup's behavior towards duplicate tags/pairs at the exact
        # same location.
        # 'auto': Calculate the maximum tags at the exact same location based on
        # binomial distribution.
        # integer value: Keep at most that much reads.
        # 'all': Keep all duplicates.
        keep_dup = 1,
        # Set verbose level.
        # 0: only show critical message.
        # 1: show additional warning message.
        # 2: show process information.
        # 3: show debug messages.
        # If you want to know where are the duplicate reads, use 3.
        # Default: 2
        verbose = 2,
    threads: 1  # Multithreading not supported.
    benchmark:
        repeat('benchmarks/macs2_filterdup/{sample}.tsv', 1)
    log: 'logs/macs2_filterdup/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/macs2/filterdup'

#
# Peak calling for each replicate
#
def macs2_callpeak_for_replicate_input(wildcards):
    # wildcards: name, i (= replicate number)
    control1, control2, treat = namerep2ct[(wildcards.name, 'rep' + wildcards.i)]

    control1 = str(result_dir / '03_macs2' / f'{control1}.sorted.filterdup.bed')
    control2 = str(result_dir / '03_macs2' / f'{control2}.sorted.filterdup.bed')
    treat = str(result_dir / '03_macs2' / f'{treat}.sorted.filterdup.bed')
    
    return {
        'control': [control1, control2],
        'treat': treat
    }
    
rule macs2_callpeak_narrow_for_replicate:
    input: unpack(macs2_callpeak_for_replicate_input)
    output:
        result_dir / '04_macs2_callpeak' / 'narrow_rep{i}' / '{name}_peaks.narrowPeak',
        result_dir / '04_macs2_callpeak' / 'narrow_rep{i}' / '{name}_treat_pileup.bdg',
        result_dir / '04_macs2_callpeak' / 'narrow_rep{i}' / '{name}_control_lambda.bdg',
    params:
        outdir = str(result_dir / '04_macs2_callpeak' / 'narrow_rep{i}')
    shadow: 'shallow'
    shell:
        'macs2 callpeak -f BEDPE --tempdir . '
        '-t {input.treat} -c {input.control} -n {wildcards.name} --outdir {params.outdir} --bdg ' # --bdg to get pileup.
        '--min-length 250'

rule macs2_callpeak_broad_for_replicate:
    input: unpack(macs2_callpeak_for_replicate_input)
    output:
        result_dir / '04_macs2_callpeak' / 'broad_rep{i}' / '{name}_peaks.broadPeak'
    params:
        outdir = str(result_dir / '04_macs2_callpeak' / 'broad_rep{i}')
    shadow: 'shallow'
    shell:
        'macs2 callpeak -f BEDPE --tempdir . '
        '--broad -t {input.treat} -c {input.control} -n {wildcards.name} --outdir {params.outdir} '
        '--min-length 250'

rule merge_narrow_broad_peaks_for_replicate:
    input:
        result_dir / '04_macs2_callpeak' / 'narrow_rep{i}' / '{name}_peaks.narrowPeak',
        result_dir / '04_macs2_callpeak' / 'broad_rep{i}' / '{name}_peaks.broadPeak',
    output:
        result_dir / '04_macs2_callpeak' / 'merged_rep{i}' / '{name}_merged_peaks.bed',
    shell:
        'cat {input} | cut -f1-3 | bedtools sort -i stdin | bedtools merge > {output}'

rule macs2_bdgcmp_ppois_for_replicate:
    input:
        treat = result_dir / '04_macs2_callpeak' / 'narrow_rep{i}' / '{name}_treat_pileup.bdg',
        control = result_dir / '04_macs2_callpeak' / 'narrow_rep{i}' / '{name}_control_lambda.bdg',
    output:
        result_dir / '05_macs2_bdg' / 'narrow_rep{i}' / '{name}_ppois.bdg'
    params:
        outdir = str(result_dir / '05_macs2_bdg' / 'narrow_rep{i}')
    shell:
        'macs2 bdgcmp -t {input.treat} -c {input.control} -o {output} -m ppois'

#
# Peak calling for pooled replicates
#
def macs2_callpeak_for_pooled_input(wildcards):
    # wildcards: name

    control1, control2, treat1 = namerep2ct[(wildcards.name, 'rep1')]
    control1, control2, treat2 = namerep2ct[(wildcards.name, 'rep2')]

    control1 = str(result_dir / '03_macs2' / f'{control1}.sorted.filterdup.bed')
    treat1 = str(result_dir / '03_macs2' / f'{treat1}.sorted.filterdup.bed')
    control2 = str(result_dir / '03_macs2' / f'{control2}.sorted.filterdup.bed')
    treat2 = str(result_dir / '03_macs2' / f'{treat2}.sorted.filterdup.bed')
    
    return {
        'control': [control1, control2],
        'treat': [treat1, treat2]
    }
    
rule macs2_callpeak_narrow_for_pooled:
    input: unpack(macs2_callpeak_for_pooled_input)
    output:
        result_dir / '04_macs2_callpeak' / 'narrow_pooled' / '{name}_peaks.narrowPeak',
        result_dir / '04_macs2_callpeak' / 'narrow_pooled' / '{name}_treat_pileup.bdg',
        result_dir / '04_macs2_callpeak' / 'narrow_pooled' / '{name}_control_lambda.bdg',
    params:
        outdir = str(result_dir / '04_macs2_callpeak' / 'narrow_pooled')
    shadow: 'shallow'
    shell:
        'macs2 callpeak -f BEDPE --tempdir . '
        '-t {input.treat} -c {input.control} -n {wildcards.name} --outdir {params.outdir} --bdg '
        '--min-length 250'

rule macs2_callpeak_broad_for_pooled:
    input: unpack(macs2_callpeak_for_pooled_input)
    output:
        result_dir / '04_macs2_callpeak' / 'broad_pooled' / '{name}_peaks.broadPeak'
    params:
        outdir = str(result_dir / '04_macs2_callpeak' / 'broad_pooled')
    shadow: 'shallow'
    shell:
        'macs2 callpeak -f BEDPE --tempdir . '
        '--broad -t {input.treat} -c {input.control} -n {wildcards.name} --outdir {params.outdir} '
        '--min-length 250'

rule merge_narrow_broad_peaks_for_pooled:
    input:
        result_dir / '04_macs2_callpeak' / 'narrow_pooled' / '{name}_peaks.narrowPeak',
        result_dir / '04_macs2_callpeak' / 'broad_pooled' / '{name}_peaks.broadPeak',
    output:
        result_dir / '04_macs2_callpeak' / 'merged_pooled' / '{name}_merged_peaks.bed',
    shell:
        'cat {input} | cut -f1-3 | bedtools sort -i stdin | bedtools merge > {output}'

rule macs2_bdgcmp_ppois_for_pooled:
    input:
        treat = result_dir / '04_macs2_callpeak' / 'narrow_pooled' / '{name}_treat_pileup.bdg',
        control = result_dir / '04_macs2_callpeak' / 'narrow_pooled' / '{name}_control_lambda.bdg',
    output:
        result_dir / '05_macs2_bdg' / 'narrow_pooled' / '{name}_ppois.bdg'
    params:
        outdir = str(result_dir / '05_macs2_bdg' / 'narrow_pooled')
    shell:
        'macs2 bdgcmp -t {input.treat} -c {input.control} -o {output} -m ppois'
