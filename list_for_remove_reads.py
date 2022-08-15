#! /usr/bin/env python

import os
import glob
import pandas as pd
import itertools

import sys
path_old_txt = sys.argv[1]
path_new_txt = sys.argv[2]

#"/Users/unnost/workwork/prespacers_task/probnaya/Results_fragseq/ids_uni.txt"

dff = pd.read_csv(path_old_txt, header=None)
#dff.head()
dff.drop_duplicates(subset =0, keep =False, inplace = True)
ids = list(dff[0])
fp = open(path_new_txt, "a")
for item in ids:
    fp.write("%s\n" % item)
fp.close()