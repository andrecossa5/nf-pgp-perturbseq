#!/usr/bin/python

import sys
import pandas as pd
import numpy as np
import pandas as pd
import dask.dataframe as dd
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt


##


# Utils
def _rev(x):
    d = {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'N':'N'}
    x = list(x)
    rev_x = []
    for i in range(len(x)):
        rev_x.append(d[x[i]])
    return ''.join(rev_x)[::-1]


##


def rev_complement(s):
    return s.map(_rev)


##


def to_numeric(X):
    return np.select([X=='A', X=='T', X=='C', X=='G'], [1,2,3,4], default=0)


##


# Args
sample = sys.argv[1]
path_bulk = sys.argv[2]
path_sc = sys.argv[3]
method = sys.argv[4]

# sample = 'AA'
# path_bulk = '/Users/IEO5505/Desktop/PD/tmp_bulk/clonal_prevalences.csv'
# path_sc = '/Users/IEO5505/Desktop/PD/tmp_sc/GBC_read_elements.tsv'


##


# Bulk prevalences
bulk = pd.read_csv(path_bulk, index_col=0)

# Read single-cell read elements, reverse-complement and read_count all sc GBCs
sc_df = dd.read_csv(path_sc, sep='\t', header=None)
sc_df.columns = ['name', 'CBC', 'UMI', 'GBC']
sc_df['GBC'] = sc_df['GBC'].map_partitions(rev_complement, meta=('GBC', 'str'))
sc = sc_df['GBC'].value_counts().compute()

##

# Map sc GBCs to bulk GBCs: if hamming <=3, correct sc GBCs to bulk GBCs
sc_A = to_numeric(np.vstack(sc.index.map(lambda x: np.array(list(x)))))
bulk_A = to_numeric(np.vstack(bulk.index.map(lambda x: np.array(list(x)))))
D = pairwise_distances(sc_A, bulk_A, metric='hamming', n_jobs=8) * sc_A.shape[1]
df_correction = (
    sc.to_frame('read_count')
    .assign(
        final_GBC=[ bulk.index[i] for i in D.argmin(axis=1) ],
        hamming=D.min(axis=1),
        to_correct=np.sum(D<=3, axis=1)>0
    )
)
test = ~df_correction['to_correct']
df_correction.loc[test, 'final_GBC'] = df_correction.index[test]
mapping = df_correction['final_GBC'].to_dict()
sc_df['GBC'] = sc_df['GBC'].map(mapping)                # Correct GBC here!

##

# Compute sc CBC-GBCs combos
df_combos = (
    sc_df.groupby(['CBC', 'GBC']).size().compute().to_frame('read_counts')
    .join(
        sc_df.groupby(['CBC', 'GBC', 'UMI']).size().compute().to_frame('umi_counts')
    )
    .reset_index()
)
df_combos['coverage'] = df_combos['read_counts'] / df_combos['umi_counts']
df_combos.to_csv('CBC_GBC_combos.tsv.gz', sep='\t')


##


# Cell assignment

# Cell assignment, only unique combos
t = 10
test = (df_combos['read_counts'] >= 30) & \
       (df_combos['umi_counts'] >= 5) & \
       (df_combos['coverage'] >= t)
df_combos['status'] = np.where(test, 'supported', 'unsupported')

##

# Combinations support plot
fig, ax = plt.subplots()
x = np.log10(df_combos['read_counts'])
y = np.log10(df_combos['umi_counts'])
ax.plot(x[df_combos['status'] == 'supported'], 
        y[df_combos['status'] == 'supported'], '.', label='assigned', color='blue')
ax.plot(x[df_combos['status'] == 'unsupported'],
        y[df_combos['status'] == 'unsupported'], '.', label='not-assigned', color='grey')
ax.set(
    title='CBC-GBC combination status', 
    xlabel='log10_read_counts', 
    ylabel='log10_umi_counts'
)
ax.legend()
fig.savefig('CBC_GBC_combo_status.png')

##


# Cell assignment and clone calling

if method == 'unique_GBC':

    # Uniquely assigned cells
    uniquely_assigned = (
        df_combos
        .groupby(['CBC', 'status'])['GBC']
        .nunique().reset_index()
        .query('status=="supported" and GBC==1')
        ['CBC'].unique()
    )
    cell_df = (
        df_combos
        .query('CBC in @uniquely_assigned and status=="supported"')
        [['CBC', 'GBC']].drop_duplicates()
        .set_index('CBC')
        .rename(columns={'GBC':'clone'})
    )
    cell_df.to_csv('cells_summary_table.csv')
    clones_df = (
        cell_df['clone']
        .value_counts(normalize=True)
        .to_frame('clonal_prevalence')
    )
    clones_df.to_csv('clones_summary_table.csv')

    ## 

elif method == 'unique_combos':

    # Unique combos
    M = (
        df_combos
        .query('status=="supported"')
        .pivot_table(index='CBC', columns='GBC', values='umi_counts')
    )
    M[M.isna()] = 0
    cells_with_unique_combos = M.apply(lambda x: frozenset(M.columns[x>0]), axis=1)

    clones_df = cells_with_unique_combos.value_counts(normalize=True)
    clones_df = (
        clones_df
        .to_frame('clonal_prevalence').reset_index()
        .rename(columns={'index':'GBC_set'})
        .assign(clone=[ f'C{i}_{sample}' for i in range(clones_df.shape[0]) ])
        .set_index('clone')
    )
    clones_df.to_csv('clones_summary_table.csv')
    cell_df = (
        cells_with_unique_combos
        .map(lambda x: clones_df.index[clones_df['GBC_set']==x][0])
        .to_frame('clone')
    )
    cell_df.to_csv('cells_summary_table.csv')

else:

    raise ValueError(
        f'{method} not supported. Supported methods are unique_combos and unique_GBC.'
    )

##