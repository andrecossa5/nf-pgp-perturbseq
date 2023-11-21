#!/usr/bin/python

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dask.dataframe as dd

##


def cell_assignment(df, t):
    """
    Assing cells to GBCs.
    """
       
    test = (df['read_counts'] > 30) & (df['umi_counts'] > 3) & (df['coverage'] > t)
    df['status'] = np.where(test, 1, 0)

    # Cell MOI tally 
    cell_MOI = (
        df.loc[:, ['CBC', 'GBC', 'status']]
        .query('status == 1')
        .groupby('CBC')
        .sum()
        .rename(columns={'status': 'multiplicity'})
    )
    # Filter uniquely assigned cells and get their unique GBC
    uniquely_assigned_cells = cell_MOI.query('multiplicity == 1').index

    return df, uniquely_assigned_cells


##


# Paths
path_elements = sys.argv[1]
path_aligned_reads = sys.argv[4]

# Read GBC-containing reads components, and reads aligned to the bulk reference
df = dd.read_csv(path_elements, sep='\t', header=None, dtype='string')
df.columns = ['read_name', 'CBC', 'GBC', 'UMI']
aligned_names = pd.read_csv(path_aligned_reads, sep='\t', header=None, dtype='string')[0]

# Filter for only reads aligned to the bulk reference, write all CBC-UMI-GBC table
df = df.loc[df['read_name'].isin(aligned_names),:]
df.compute().to_csv('CBC_UMI_GBC_by_read.tsv.gz', sep='\t', index=False)

# Compute combos table
df_combos = (
    df.groupby(['CBC', 'GBC']).size().compute().to_frame('read_counts')
    .join(
        df.groupby(['CBC', 'GBC', 'UMI']).size().compute().to_frame('umi_counts')
    )
    .reset_index()
)
df_combos['coverage'] = df_combos['read_counts'] / df_combos['umi_counts']
df_combos.to_csv('CBC_GBC_combos.tsv.gz', sep='\t')


##


# Cell classification: from Adamson et al. 2016
L = []
tresholds = np.arange(5, 100, 5)
for t in tresholds:
    _, assigned_cells = cell_assignment(df_combos, t)
    L.append(len(assigned_cells))
optimal_t = tresholds[np.argmax(L)]
df_combos, uniquely_assigned_cells = cell_assignment(df_combos, optimal_t)

# Cell assignment plot
fig, ax = plt.subplots()
x = np.log10(df_combos['read_counts'])
y = np.log10(df_combos['umi_counts'])
ax.plot(x[df_combos['status'] == 1], y[df_combos['status'] == 1], '.', label='assigned', color='blue')
ax.plot(x[df_combos['status'] == 0], y[df_combos['status'] == 0], '.', label='not-assigned', color='grey')
ax.set(title='CBC-GBC combination status', xlabel='log10_read_counts', ylabel='log10_umi_counts')
ax.legend()
fig.savefig('CBC_GBC_combo_status.png')

# Compute summary tables: cells 
df_cells = (
    df_combos
    .query('CBC in @uniquely_assigned_cells and status == 1')
    .drop(columns='status')
    .set_index('CBC')
)

# Assert we have taken the correct ones and save
assert df_cells.index.size == uniquely_assigned_cells.size
df_cells.to_csv('cells_summary_table.csv')

# Compute summary tables: clones
df_clones = df_cells.reset_index().groupby('GBC').size().to_frame('n_cells')
df_clones['prevalence'] = df_clones['n_cells'] / df_clones['n_cells'].sum()
df_clones = df_clones.sort_values('prevalence', ascending=False)
df_clones.to_csv('clones_summary_table.csv')


##
