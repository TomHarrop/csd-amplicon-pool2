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

logging.info(f'Reading {my_bam}')
bamfile = pysam.AlignmentFile(my_bam, 'rb')
logging.info(f'{my_bam} contains {bamfile.count()} reads')

logging.info(f'Finding reads mapped to {hvr_chr}:{hvr_start}')
start_iter = bamfile.fetch(hvr_chr, hvr_start)

logging.info(f'Finding reads mapped to {hvr_chr}:{hvr_stop}')
stop_iter = bamfile.fetch(hvr_chr, hvr_stop)
stop_ids = [x.qname for x in stop_iter]
logging.info(f'{len(stop_ids)} reads mapped to {hvr_chr}:{hvr_stop}')

logging.info(f'Finding reads spanning {hvr_chr}:{hvr_start}-{hvr_stop}')
span_ids = sorted(set(x.qname for x in start_iter
                      if (x.qname in stop_ids
                          and len(x.seq) < 500
                          and x.mapq > 30)))

logging.info(f'Writing {len(span_ids)} ids to {id_file}')
with open(id_file, 'wt') as f:
    f.write('\n'.join(span_ids))
    f.write('\n')
