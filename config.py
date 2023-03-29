#!/usr/bin/env python

import os
from filo import pfam, analysis

# defining variables

dicio = {'specie': ['code', 'abreviation']}
pfam_path = ''
regex = ''
seq_path = ''
seq_list = os.listdir(seq_path)

print('\n Retriving database')
pfam(seq_list=seq_list, seq_path=seq_path, dicio=dicio, pfam_path=pfam_path, regex=regex)

print('\n Starting gene prospecting analyses')
outdir = ''
domain = []
analysis(outdir=outdir, domain=domain, seq_list=seq_list, seq_path=seq_path, dicio=dicio, regex=regex)
