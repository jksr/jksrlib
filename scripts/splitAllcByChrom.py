#!/usr/bin/env python

import sys
import subprocess
import os

allcin = sys.argv[1]
outprefix = sys.argv[2]
chrs = sys.argv[3] if len(sys.argv)==4 else None

if chrs is None:
	chrs =[ x.decode("utf-8")  for x in  subprocess.check_output(['tabix','-l',allcin]).split() ]

	tmp = []
	for chrom in chrs:
		which = chrom.replace('chr','')
		if which.isnumeric() or which in ['X','Y','M']:
			tmp.append(chrom)
	chrs = tmp


for chrom in chrs:
	outallc = f'{outprefix}.{chrom}.allc.tsv.gz'
	os.system(f'tabix {allcin} {chrom} | bgzip  > {outallc}')
	os.system(f'tabix -f -b 2 -e 2 -s 1 {outallc}')


