#!/usr/bin/python

import sys
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d 
import matplotlib.pyplot as plt


##


# Helper functions
def rev_complement(x):
    d = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}
    x = list(x)
    rev_x = []
    for i in range(len(x)):
        rev_x.append(d[x[i]])
    return ''.join(rev_x)[::-1]

##


def fit_trend(df):
    """
    Plot the Log f, cell 
    """
    x = df['log_10_read_count']
    y = df['log_10_n_cells']
    f = interp1d(x, y)

    fig, ax = plt.subplots()
    ax.plot(x, y, 'o')
    ax.plot([x[0], x[-1]], [f(x)[0], f(x)[-1]], '--')
    ax.set(xlabel='Log10_read_count', ylabel='Log10_n_cells', title='Spikeins fit')
    
    return f, fig


##


# Paths 
path_spikeins_table = sys.argv[1]
path_counts = sys.argv[2]

# Read tables
df_spike = pd.read_csv(path_spikeins_table, index_col=0)
df_spike.index = df_spike.index.map(lambda x: rev_complement(x))
df_counts = pd.read_csv(path_counts, sep='\t', index_col=0)

assert np.sum(df_spike.index.isin(df_counts.index)) >= 2 # At least two GBCs are needed to make the interpolation. Better all of them though...

# Fit spikeins trend
spikeins = df_spike.index[df_spike.index.isin(df_counts.index)].to_list()
df_spike = df_spike.loc[spikeins, :].join(df_counts.loc[spikeins, :])
df_spike['log_10_n_cells'] = np.log10(df_spike['n_cells'])
df_spike['log_10_read_count'] = np.log10(df_spike['read_count'])
f, fig = fit_trend(df_spike)
# Save spikins fit
fig.savefig('spikeins_fit.png')

# Remove spikeins from putative GBCs 
is_spike = df_counts.index.isin(df_spike.index)
df_counts = df_counts[~is_spike]

# Remove GBCs falling outside the spikeins (log10) interpolation range
df_counts['log_10_read_count'] = np.log10(df_counts['read_count'])
min_ = df_spike['log_10_read_count'][-1]
max_ = df_spike['log_10_read_count'][0]
is_in_range = (df_counts['log_10_read_count'] >= min_) & (df_counts['log_10_read_count'] <= max_)
df_counts = df_counts[is_in_range]

# Fit trend on their log10 read_count values, returning log10 n_cells
df_counts['log_10_n_cells'] = f(df_counts['log_10_read_count'])
# n_cells
df_counts['n_cells'] = np.round(10 ** df_counts['log_10_n_cells'])
df_counts['cellular_prevalence'] = df_counts['n_cells'] / df_counts['n_cells'].sum()

# Save
df_counts = df_counts.loc[:, ~df_counts.columns.str.contains('log_10')]
df_counts.to_csv('clonal_prevalences.csv')
df_counts.index.to_series().to_csv('good_GBCs_bulk.txt', index=False, header=False)


##