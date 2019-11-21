#!/usr/bin/env python3

import pysam
import logging

# set up log
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    filename=snakemake.log[0],
    level=logging.INFO)


my_bam = snakemake.input['bam']
# my_bam = 'output/020_mapped/minion/WS20_81_sorted.bam'
id_file = snakemake.output['ids']
hvr_chr = 'NC_037640.1'
hvr_start = 11771976
hvr_stop = 11772216
hvr_size = hvr_stop - hvr_start

logging.info(f'Reading {my_bam}')
bamfile = pysam.AlignmentFile(my_bam, 'rb')
logging.info(f'{my_bam} contains {bamfile.count()} reads')

logging.info(f'Finding reads mapped to {hvr_chr}:{hvr_start}-{hvr_stop}')
span_iter = bamfile.fetch(hvr_chr, hvr_start, hvr_stop)

logging.info(f'Extracting reads with > 90% coverage of HVR exon')
kept_ids = sorted(set(x.qname for x in span_iter
                      if (x.get_overlap(hvr_start, hvr_stop) / hvr_size > 0.9
                          and len(x.seq) < 500
                          and len(x.seq) > 370
                          and x.mapq > 30)))
logging.info(f'{len(kept_ids)} reads kept')

logging.info(f'Writing ids to {id_file}')
with open(id_file, 'wt') as f:
    f.write('\n'.join(kept_ids))
    f.write('\n')
