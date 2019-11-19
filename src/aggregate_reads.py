#!/usr/bin/env python3

import fileinput
import logging

# set up log
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    filename=snakemake.log[0],
    level=logging.DEBUG)

# catch files from snakemake
input_list = snakemake.input
fq = snakemake.output['fq']

logging.debug('input_list')
logging.debug(f'{input_list}')
logging.debug('fq')
logging.debug(f'{fq}')

with open(fq, 'wt') as f:
    for line in fileinput.input(input_list):
        f.write(line)
