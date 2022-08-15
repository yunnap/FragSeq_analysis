#! /usr/bin/env python

import os
import glob
import pandas as pd
import itertools

import sys
in_bam = sys.argv[1]
out_bam = sys.argv[2]

import pysam

inbam = pysam.AlignmentFile(in_bam)
outbam = pysam.AlignmentFile(out_bam, "wb", template=inbam)

for read in inbam.fetch(until_eof=True):
    if 'I' not in read.cigarstring and 'D' not in read.cigarstring:
        outbam.write(read)

        