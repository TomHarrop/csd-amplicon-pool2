#!/usr/bin/env python3

import csv
import multiprocessing
import pathlib
import re

#############
# FUNCTIONS #
#############

def aggregate_raw_reads(wildcards):
    my_bc = indiv_to_bc[wildcards.indiv]
    my_bc_str = re.findall('\d+', my_bc)[0]
    my_bc_folder = 'barcode' + my_bc_str
    my_path = f'data/reads/{wildcards.run}/pass/{my_bc_folder}'
    my_files = pathlib.Path(my_path).glob('*.fastq')
    return([str(x) for x in my_files])

###########
# GLOBALS #
###########

sample_key = 'data/Pool 2 barcodes.csv'

bioconductor_container = 'shub://TomHarrop/singularity-containers:bioconductor_3.9'
biopython_container = 'shub://TomHarrop/singularity-containers:biopython_1.73'
freebayes_container = 'shub://TomHarrop/singularity-containers:freebayes_1.2.0'
minimap_container = 'shub://TomHarrop/singularity-containers:minimap2_2.11r797'
sambamba_container = 'shub://TomHarrop/singularity-containers:sambamba_0.6.9'
samtools_container = 'shub://TomHarrop/singularity-containers:samtools_1.9'

########
# MAIN #
########

# map indiv to barcode
indiv_to_bc = {}
with open(sample_key, 'rt') as f:
    csv_reader = csv.reader(f)
    next(csv_reader)
    for row in csv_reader:
        my_bc = row[1]
        my_bc_str = re.findall('\d+', my_bc)[0]
        my_id = row[2] + '_' + my_bc_str
        indiv_to_bc[my_id] = my_bc

bc_to_indiv = {indiv_to_bc[x]: x for x in indiv_to_bc.keys()}

# exclude barcode 25
all_indivs = [x for x in indiv_to_bc.keys() if indiv_to_bc[x] != 'BC25']

all_bc = sorted(set(glob_wildcards('data/reads/minion/pass/{bc}/{id}.fastq').bc))

#########
# RULES #
#########


rule target:
    input:
        expand('output/050_derived-alleles/{run}/all-indivs_aa.fa',
               run=['flongle']),     # not enough RAM to basecall minion run
        'output/020_mapped/flongle/merged.bam'


rule translate_consensus:
    input:
        'output/050_derived-alleles/{run}/consensus_all-indivs.fa'
    output:
        'output/050_derived-alleles/{run}/all-indivs_aa.fa'
    singularity:
        biopython_container
    script:
        'src/translate_consensus.py'

rule combine_cds:
    input:
        expand('output/050_derived-alleles/{{run}}/{indiv}_consensus_condensed.fa',
               indiv=all_indivs)
    output:
        'output/050_derived-alleles/{run}/consensus_all-indivs.fa'
    singularity:
        samtools_container
    shell:
        'cat {input} > {output}'

rule condense_cds:
    input:
        'output/050_derived-alleles/{run}/{indiv}_consensus.fa'
    output:
        temp('output/050_derived-alleles/{run}/{indiv}_consensus_condensed.fa')
    params:
        header = '>{indiv}'
    singularity:
        samtools_container
    shell:
        'echo "{params.header}" > {output} ; '
        'grep -v "^>" {input} >> {output}'


rule extract_derived_cds:
    input:
        fa = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna',
        regions = 'output/040_variant-annotation/csd_regions.txt',
        vcf = 'output/040_variant-annotation/{run}_csd_reheadered.vcf.gz'
    output:
        temp('output/050_derived-alleles/{run}/{indiv}_consensus.fa')
    log:
        'output/logs/050_derived-alleles/{run}_{indiv}_consensus.log'
    singularity:
        samtools_container
    shell:
        'samtools faidx '
        '{input.fa} '
        '$(cat {input.regions}) '
        '2> {log} '
        '| '
        'bcftools consensus '
        '-s {wildcards.run}_{wildcards.indiv} '
        '-H 1 '
        '{input.vcf} '
        '> {output} '
        '2>> {log}'


rule filter_csd_variants:
    input:
        gff = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.gff',
        vcf = 'output/030_freebayes/{run}/variants.vcf.gz',
        tbi = 'output/030_freebayes/{run}/variants.vcf.gz.tbi',
        fa = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna',
    output:
        coding = 'output/040_variant-annotation/{run}_coding.Rds',
        csd = 'output/040_variant-annotation/{run}_csd.vcf',
    log:
        'output/logs/040_variant-annotation/{run}_filter_csd_variants.log'
    singularity:
        bioconductor_container
    script:
        'src/filter_csd_variants.R'


rule freebayes:
    input:
        bam = expand('output/020_mapped/{{run}}/{indiv}_sorted.bam',
                     indiv=all_indivs),
        fa = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'
    output:
        vcf = 'output/030_freebayes/{run}/variants.vcf'
    params:
        region = 'NC_037640.1:11771679-11781139'
    log:
        'output/logs/030_freebayes/{run}-freebayes.log'
    singularity:
        freebayes_container
    shell:
        'freebayes '
        '--region {params.region} '
        '-f {input.fa} '
        '{input.bam} '
        '> {output} '
        '2> {log}'


rule merge_bam: # for visualisation
    input:
        expand('output/020_mapped/{{run}}/{indiv}.sam',
               indiv=all_indivs)
    output:
        'output/020_mapped/{run}/merged.bam'
    log:
        'output/logs/020_mapped/{run}/merge.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        sambamba_container
    shell:
        'sambamba merge '
        '--nthreads={threads} '
        '--compression-level=9 '
        '{output} '
        '{input} '
        '2> {log} '
        '; '
        'sambamba index {output} 2>> {log}'

rule sort_sam:
    input:
        'output/020_mapped/{run}/{indiv}.sam'
    output:
        bam = 'output/020_mapped/{run}/{indiv}_sorted.bam'
    log:
        'output/logs/020_mapped/{run}/{indiv}_sort.log'
    threads:
        1
    singularity:
        sambamba_container
    shell:
        'sambamba view '
        '{input} '
        '-f "bam" '
        '--sam-input '
        '-l 0 '
        '2> {log} '
        '| '
        'sambamba sort '
        '-o {output.bam} '
        '-l 9 '
        '-t {threads} '
        '/dev/stdin '
        '2>> {log} '
        '; '
        'sambamba index '
        '{output.bam} '
        '2>> {log}'

rule map_to_genome:
    input:
        fq = 'output/010_raw/{run}/{indiv}.fq',
        ref = 'output/010_raw/honeybee_ref.mmi'
    output:
        temp('output/020_mapped/{run}/{indiv}.sam')
    params:
        rg = '\'@RG\\tID:{run}_{indiv}\\tSM:{run}_{indiv}\''
    log:
        'output/logs/020_mapped/{run}/{indiv}.log'
    threads:
        1
    singularity:
        minimap_container
    shell:
        'minimap2 '
        '-t {threads} '
        '-ax map-ont '
        '-N 1 '
        '-R {params.rg} '
        '{input.ref} '
        '{input.fq} '
        '> {output} '
        '2> {log}'

rule aggregate_reads:
    input:
        aggregate_raw_reads
    output:
        'output/010_raw/{run}/{indiv}.fq'
    shell:
        'cat {input} > {output}'

rule prepare_ref:
    input:
        'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'
    output:
        'output/010_raw/honeybee_ref.mmi'
    log:
        'output/logs/010_raw/prepare_ref.log'
    priority:
        100
    threads:
        3
    singularity:
        minimap_container
    shell:
        'minimap2 '
        '-x map-ont '
        '-d {output} '
        '{input} '
        '2> {log}'


# generic csd rules
rule extract_csd_regions:
    input:
        gff = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.gff'
    output:
        regions = 'output/040_variant-annotation/csd_regions.txt'
    log:
        'output/logs/040_variant-annotation/extract_csd_regions.log'
    singularity:
        bioconductor_container
    script:
        'src/extract_csd_regions.R'


# generic index rule
rule index_vcf:
    input:
        'output/{folder}/{file}.vcf'
    output:
        gz = 'output/{folder}/{file}.vcf.gz',
        tbi = 'output/{folder}/{file}.vcf.gz.tbi'
    log:
        'output/logs/{folder}/{file}_index-vcf.log'
    singularity:
        samtools_container
    shell:
        'bgzip -c {input} > {output.gz} 2> {log} '
        '; '
        'tabix -p vcf {output.gz} 2>> {log}'

# generic reheader rule
rule reheader_vcf:
    input:
        'output/{folder}/{file}.vcf'
    output:
        'output/{folder}/{file}_reheadered.vcf'
    singularity:
        samtools_container
    shell:
        'grep -v "^##FILTER=All filters passed" {input} > {output}'


