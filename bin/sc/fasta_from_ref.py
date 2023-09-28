#!/usr/bin/python

import os
import sys
import pandas as pd


##


# Paths
path_sample_folder = sys.argv[1]

# Read clone prevalences
df = pd.read_csv(
    os.path.join(path_sample_folder, 'clonal_prevalences.csv'), 
    index_col=0
)
GBCs = df.index.to_list()

# Write as fasta
with open('GBC_reference.fa', 'w') as f:
    for i, x in GBCs:
        f.write(f'>seq{i}\n')
        f.write(f'{x}\n')


##
