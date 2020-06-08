#! /usr/bin/env python3

import os
import glob
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--count', metavar='reads_count_folder', dest='c',
                    type=str, required=True)
parser.add_argument('-l', '--length', metavar='gene_lengths', dest='l',
                    type=str, required=True)
parser.add_argument('-o', '--output', metavar='output_RPKM_table', dest='o',
                    type=str, required=True)
args = parser.parse_args()


# parse gene lengths
dict_l = {}
with open(args.l, 'r') as fi:
        for line in fi:
            (key, val) = line.strip().split('\t')
            dict_l[key] = val

# parse *.count files
for f in glob.glob(os.path.join(args.c ,'*.count')):
    pd.read_csv(f,sep='\t')

df_rcount = pd.concat([pd.read_csv(f,sep='\t',index_col=0, names = ['gene', os.path.split(f)[1]]) for f in glob.glob(os.path.join(args.c ,'*.count'))], axis=1)
df_glen = pd.read_csv(args.l, sep='\t', index_col=0, names=['gene', 'length'])
df_cat = pd.concat([df_rcount, df_glen], axis=1, join='inner')
for index, row in df_cat.iteritems():
    if index != 'length':
        nr = df_cat[index].sum()
        df_cat[index] = nr / df_cat['length'] / df_cat[index] * 1e9
df_cat = df_cat.drop(columns='length')
df_cat.replace(np.inf, 0, inplace=True)
df_cat.to_csv(args.o, sep='\t') 


