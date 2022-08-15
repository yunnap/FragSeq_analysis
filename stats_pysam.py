#! /usr/bin/env python

import os
import glob
import pandas as pd
import itertools

import sys
path_bam = sys.argv[1]
path_stat_txt = sys.argv[2]

import pysam

f_bam = pysam.AlignmentFile(path_bam)
stat = f_bam.get_index_statistics()
#stat = list(stat)



with open(path_stat_txt, "a") as file:
    #print(stat, file=file, sep="\n")
    print("\n".join(map(str, stat)), file=file)



oligs = ['20','30','31','32','33','34','35','36','37','38','39','40','50','60']

try:
	sum_oligs_map = 0
	for elem in stat:
	    el=str(elem).split('(')[1].split(',')[0].split('=')[1].split("'")[1]
	    if el in oligs:
	        num_map = int(str(elem).split('(')[1].split(',')[1].split('=')[1].split("'")[0])
	        sum_oligs_map += num_map

	with open(path_stat_txt, "a") as file:
	    print("", file=file)
	    print("oligs_mapped =", sum_oligs_map, file=file)
	        	

except:
	with open(path_stat_txt, "a") as file:
		print("", file=file)
		print('cannot count the sum of the reads mapped on oligs', file=file)

#fp = open(path_stat_txt, "a")
#for item in stat:
#    fp.write("%s\n" % item)
#fp.close()


