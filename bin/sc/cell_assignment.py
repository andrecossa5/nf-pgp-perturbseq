#!/usr/bin/python

# Cell assignment script

########################################################################

# Parsing CLI args 
 
# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='cell_assignment',
    description=
    """
    Script for clone calling and cell assignment.
    """
)


# Input
my_parser.add_argument(
    '--sample', 
    type=str,
    default=None,
    help='Sample name. Default: None.'
)

# Output
my_parser.add_argument(
    '--path_bulk', 
    type=str,
    default=None,
    help='Path to input bulk reference. Default: None.'
)

# treshold
my_parser.add_argument(
    '--path_sc',
    type=str,
    default=None,
    help='Path to input sc GBC reads elements. Default: None.'
)

# treshold
my_parser.add_argument(
    '--method',
    type=str,
    default='unique_combos',
    help='Method for clone calling and cell assignment. Default: unique_combos.'
)

# treshold
my_parser.add_argument(
    '--ncores',
    type=int,
    default=8,
    help='n cores for pairwise distances calculation. Default: 8.'
)

# Spikeins
my_parser.add_argument(
    '--bulk_sc_treshold',
    type=str,
    default=1,
    help='''
    Hamming distance treshold to consider a sc GBC a "degenerate" sequence with respect 
    to a bulk reference sequence. Default: 1.
    '''
)

# Spikeins
my_parser.add_argument(
    '--umi_treshold',
    type=str,
    default=5,
    help='Min number of UMIs to consider a CB-GBC combination supported. Default: 5.'
)

# Spikeins
my_parser.add_argument(
    '--read_treshold',
    type=str,
    default=15,
    help='Min number of reads to consider a CB-GBC combination supported. Default: 15.'
)

# Spikeins
my_parser.add_argument(
    '--coverage_treshold',
    type=str,
    default=3,
    help='Min coverage (nUMIs / nreads) to consider a CB-GBC combination supported. Default: 3.'
)


##


# Parse arguments
args = my_parser.parse_args()
sample = args.sample
path_bulk = args.path_bulk
path_sc = args.path_sc
method = args.method
ncores = args.ncores
bulk_sc_treshold = args.bulk_sc_treshold
umi_treshold = args.umi_treshold
read_treshold = args.read_treshold
coverage_treshold = args.coverage_treshold


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


# Import code
import pandas as pd
import numpy as np
import pandas as pd
import dask.dataframe as dd
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt


##


########################################################################

def main():

    # Bulk reference
    bulk = pd.read_csv(
        os.path.join(path_bulk, 'summary', 'bulk_GBC_reference.csv'),
        index_col=0, header=None
    )

    # Read single-cell read elements, reverse-complement GBCs and count reads
    sc_df = dd.read_csv(path_sc, sep='\t', header=None)
    sc_df.columns = ['name', 'CBC', 'UMI', 'GBC']
    sc_df['GBC'] = sc_df['GBC'].map_partitions(rev_complement, meta=('GBC', 'str'))
    sc = sc_df['GBC'].value_counts().compute()

    
    ##

    
    # Map sc GBCs to bulk GBCs: if hamming <=bulk_sc_treshold, correct sc GBCs to bulk GBCs
    # else, discard the read supporting that CBC-GBC combination
    sc_A = to_numeric(np.vstack(sc.index.map(lambda x: np.array(list(x)))))
    bulk_A = to_numeric(np.vstack(bulk.index.map(lambda x: np.array(list(x)))))
    D = pairwise_distances(sc_A, bulk_A, metric='hamming', n_jobs=int(ncores)) * sc_A.shape[1]
    d_correction = sc.to_frame().reset_index(drop=True)
    d_correction.columns = ['read_count']

    d_correction = (
        d_correction
        .assign(
            correct_GBC=[ bulk.index[i] for i in D.argmin(axis=1) ],
            hamming=D.min(axis=1),
        )
        .query('hamming<=@bulk_sc_treshold')
    )

    print(type(d_correction))
    print(d_correction.head())

    import sys
    sys.exit()

    # ['correct_GBC'].to_dict()
    # )

    sc_df['GBC'] = sc_df['GBC'].map(
        lambda x: d_correction[x] if x in d_correction else 'not_found'     # Correct
    )
    sc_df = sc_df.loc[lambda x: x['GBC']!='not_found']                      # Discard

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
    test = (df_combos['read_counts'] >= read_treshold) & \
           (df_combos['umi_counts'] >= umi_treshold) & \
           (df_combos['coverage'] >= coverage_treshold) 
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

########################################################################
    
# Run
if __name__ == '__main__':
    main()