#!/usr/bin/env python3

import csv
import fileinput
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

bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
bioconductor_container = 'shub://TomHarrop/singularity-containers:bioconductor_3.9'
biopython_container = 'shub://TomHarrop/singularity-containers:biopython_1.73'
clustalo = 'shub://TomHarrop/singularity-containers:clustalo_1.2.4'
freebayes_container = 'shub://TomHarrop/singularity-containers:freebayes_1.2.0'
medaka = ('shub://TomHarrop/ont-containers:medaka_7574cf1'
          '@0a7f4c632e78b4c99afbe5b52a6838d01a0cabf3')
minimap_container = 'shub://TomHarrop/singularity-containers:minimap2_2.11r797'
sambamba_container = 'shub://TomHarrop/singularity-containers:sambamba_0.6.9'
samtools_container = 'shub://TomHarrop/singularity-containers:samtools_1.9'
seqtk = ('shub://TomHarrop/seq-utils:seqtk_1.3r106'
         '@cd07278aabfe76ff52c524190ccbc5a5153a9061')


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
all_indivs = [x for x in all_indivs if x.startswith('BB')] # bb only
# all_indivs = ['BB44_60', 'WS20_81', 'TY12_28']
# all_indivs = ['BB44_60']

# clean the bamfile?
# @SQ SN:NW_020555893.1   LN:21390
# samtools view -h -F 2308 -s 0.01 -O BAM merged.bam > filtered.bam
# samtools index filtered.bam



#########
# RULES #
#########


rule target:
    input:
        expand('output/050_derived-alleles/{run}/{file}_aa.faa',
               run=[
                   # 'flongle',
                   'minion'
                    ],
               file=[
                     'all-indivs'
                     # 'drones'
                    ])

# extract and analyse results
rule align_consensus:
    input:
        'output/050_derived-alleles/{run}/{file}_aa.fa'
    output:
        aln = 'output/050_derived-alleles/{run}/{file}_aa.faa',
        dist = 'output/050_derived-alleles/{run}/{file}_aa.dist'
    log:
        'output/logs/050_derived-alleles/{run}-{file}-clustalo.fa'
    threads:
        multiprocessing.cpu_count()
    singularity:
        clustalo
    shell:
        'clustalo '
        '-i {input} '
        '--threads {threads} '
        '--dealign '
        '--full '
        '--out {output.aln} '
        '--distmat-out {output.dist} '
        '&> {log}'

rule drones_only:
    input:
        'output/050_derived-alleles/{run}/all-indivs_aa.fa'
    output:
        'output/050_derived-alleles/{run}/drones_aa.fa'
    log:
        'output/logs/050_derived-alleles/{run}-filterbyname.log'
    threads:
        1
    singularity:
        bbduk_container
    shell:
        'filterbyname.sh '
        'in={input} '
        'out={output} '
        'names=DR '
        'substring=name '
        'include=t '
        'ignorejunk=t '
        '2> {log}'

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
        expand('output/050_derived-alleles/{{run}}/{indiv}/cds_h{h}.fa',
               indiv=all_indivs,
               h=[1, 2])
    output:
        'output/050_derived-alleles/{run}/consensus_all-indivs.fa'
    singularity:
        samtools_container
    shell:
        'cat {input} > {output}'

rule condense_cds:
    input:
        'output/000_tmp/{run}/{indiv}/cds_h{h}.fa'
    output:
        'output/050_derived-alleles/{run}/{indiv}/cds_h{h}.fa'
    params:
        header = '>{indiv}_h{h}'
    singularity:
        samtools_container
    shell:
        'echo "{params.header}" > {output} ; '
        'grep -v "^>" {input} >> {output}'


rule extract_derived_cds:
    input:
        fa = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna',
        regions = 'output/005_ref/hvr_dt.txt',
        vcf = 'output/040_variant-annotation/{run}/{indiv}/filtered.vcf.gz'
    output:
        h1 = 'output/000_tmp/{run}/{indiv}/cds_h1.fa',
        h2 = 'output/000_tmp/{run}/{indiv}/cds_h2.fa'
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
        '-s {wildcards.indiv} '
        '-H 1 '
        '{input.vcf} '
        '> {output.h1} '
        '2>> {log} '
        '; '
        'samtools faidx '
        '{input.fa} '
        '$(cat {input.regions}) '
        '2>> {log} '
        '| '
        'bcftools consensus '
        '-s {wildcards.indiv} '
        '-H 1 '
        '{input.vcf} '
        '> {output.h2} '
        '2>> {log} '


# use the filtered list to extract SNPs from the non-broken VCF
rule merge_filtered_variants:  # DOESN'T WORK, RUN PER INDIV
    input:
        expand('output/040_variant-annotation/{{run}}/{indiv}/filtered.vcf.gz',
               indiv=all_indivs)
    output:
        vcf = 'output/040_variant-annotation/{run}/all_indivs.vcf',
    log:
        'output/logs/040_variant-annotation/{run}_merge.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        samtools_container
    shell:
        'bcftools merge '
        '--missing-to-ref '
        '--merge all '
        '--threads {threads} '
        '--output-type u '
        '{input} '
        '2> {log} '
        '| '
        'bcftools sort '
        '--output-file {output.vcf} '
        '--output-type v '
        '- '
        '2>> {log}'

# CONTIGL=$(awk '{printf("##contig=<ID=%s,length=%d>\\n",$1,$2);}' ../data/GCF_003254395.2_Amel_HAv3.1_genomic.fna.fai)

# grep -v "^##contig" ../output/040_variant-annotation/flongle/DR02_04/filtered.vcf |
# awk -v f="${CONTIGL}" \
#  '/^#CHROM/ { printf(f);} {print;}'  > DR02_04.vcf

# grep -v "^##contig" ../output/040_variant-annotation/flongle/DR02_94/filtered.vcf |
# awk -v f="${CONTIGL}" \
#  '/^#CHROM/ { printf(f);} {print;}'  > DR02_94.vcf


# PicardCommandLine MergeVcfs \
#     I=DR02_04.vcf \
#     I=DR02_94.vcf \
#     O=merged.vcf


rule extract_filtered_variants:
    input:
        vcf = 'output/000_tmp/{run}/{indiv}/withid.vcf.gz',
        keep = 'output/000_tmp/{run}/{indiv}/snps-to-keep.txt'
    output:
        'output/040_variant-annotation/{run}/{indiv}/filtered.vcf'
    singularity:
        samtools_container
    shell:
        'bcftools view '
        '--include ID=@{input.keep} '
        '{input.vcf} '
        '> {output}'


rule list_filtered_variants:
    input:
        'output/000_tmp/{run}/{indiv}/csd.vcf'
    output:
        'output/000_tmp/{run}/{indiv}/snps-to-keep.txt'
    singularity:
        samtools_container
    shell:
        'grep -v "^#" {input} '
        '| '
        'cut -f3 '
        '> {output}'

rule filter_csd_variants:
    input:
        txdb = 'output/005_ref/txdb.sqlite',
        fa = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna',
        vcf = 'output/000_tmp/{run}/{indiv}/withid.vcf.gz',
        tbi = 'output/000_tmp/{run}/{indiv}/withid.vcf.gz.tbi',
    output:
        csd = 'output/000_tmp/{run}/{indiv}/csd.vcf'
    log:
        'output/logs/040_variant-annotation/{run}-{indiv}-filter_csd_variants.log'
    singularity:
        bioconductor_container
    script:
        'src/filter_csd_variants.R'


rule add_snp_ids:   # otherwise annotate_variants makes them up
    input:
        'output/000_tmp/{run}/{indiv}/withsample.vcf'
    output:
        'output/000_tmp/{run}/{indiv}/withid.vcf'
    singularity:
        samtools_container
    shell:
        'bcftools annotate '
        '-I \'%CHROM\:%POS\_%REF\/%FIRST_ALT\' '
        '{input} > {output}'

# calling
rule add_sample_to_medaka:
    input:
        'output/035_medaka/{run}/{indiv}/round_1_phased.vcf'
    output:
        'output/000_tmp/{run}/{indiv}/withsample.vcf'
    singularity:
        samtools_container
    shell:
        'sed -e \'/^#CHROM/s/SAMPLE/{wildcards.indiv}/\' '
        '{input} > {output}'

rule medaka:
    input:
        bam = 'output/020_mapped/{run}/{indiv}_sorted.bam',
        fa = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'
    output:
        'output/035_medaka/{run}/{indiv}/round_1_phased.vcf'
    log:
        'output/logs/035_medaka/{run}_{indiv}.log'
    params:
        wd = 'output/035_medaka/{run}/{indiv}',
        snp_model = 'r941_min_diploid_snp',
        var_model = 'r941_min_high'
    threads:
        1
    singularity:
        medaka
    shell:
        'export TF_FORCE_GPU_ALLOW_GROWTH=true '
        '; '
        'medaka_variant '
        '-i {input.bam} '
        '-f {input.fa} '
        '-o {params.wd} '
        '-s {params.snp_model} '
        '-m {params.var_model} '
        '-t {threads} '
        '&> {log}'

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


# mapping
rule merge_bam: # for visualisation
    input:
        expand('output/020_mapped/{{run}}/{indiv}_sorted.bam',
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


rule filter_weird_reads:
    input:
        'output/010_raw/{run}/{indiv}.fq'
    output:
        pipe('{run}-{indiv}.fq')
    singularity:
        seqtk
    shell:
        'seqtk seq -C {input} >> {output}'

rule map_to_genome:
    input:
        fq = '{run}-{indiv}.fq',
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

# processing
rule generate_txdb:
    input:
        gff = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.gff',
    output:
        txdb = 'output/005_ref/txdb.sqlite'
    log:
        'output/logs/005_ref/generate_txdb.log'
    singularity:
        bioconductor_container
    script:
        'src/generate_txdb.R'

rule aggregate_reads:
    input:
        aggregate_raw_reads
    output:
        fq = 'output/010_raw/{run}/{indiv}.fq'
    log:
        'output/logs/010_raw/aggregate_reads_{run}-{indiv}.log'
    singularity:
        biopython_container
    script:
        'src/aggregate_reads.py'

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
rule extract_hvr_exon:
    input:
        gff = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.gff'
    output:
        regions = 'output/005_ref/hvr_dt.txt'
    log:
        'output/logs/005_ref/extract_hvr_exon.log'
    singularity:
        bioconductor_container
    script:
        'src/extract_hvr_exon.R'


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
