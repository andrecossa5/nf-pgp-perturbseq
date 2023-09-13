#!/usr/bin/python

# Infer clone prevalences script

########################################################################

# Parsing CLI args 

# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='infer_clone_prevalences',
    description=
    """
    Script for final inference of clone prevalence from each sample bulk GBC counts.
    """
)

# Add arguments

# Input
my_parser.add_argument(
    '-i', 
    '--input', 
    type=str,
    default=None,
    help='The path to GBC_counts.csv file. Default: . .'
)

# Output
my_parser.add_argument(
    '-o', 
    '--output', 
    type=str,
    default=os.getcwd(),
    help='Output folder path. Default: . .'
)

# Path spikeins table
my_parser.add_argument( 
    '--path_spikeins_table', 
    type=str,
    default=None,
    help='''
    Path to spikeins table. Default: None. If different from None, the script will try to use spikeins for 
    clonal inference.
    '''
)

# Cumsum
my_parser.add_argument( 
    '--cumsum', 
    type=int,
    default=75,
    help='Treshold to retain only barcodes accounting up to the <cumsum> percentile of the total GBC read_count.'
)

# Path spikeins table
my_parser.add_argument( 
    '--n_reads', 
    type=int,
    default=1000,
    help='Min n_reads to retain a GBC.'
)

# Include degenerated
my_parser.add_argument(
    '--with_degenerated', 
    action='store_true',
    help='''
        Include degenerated barcodes that were whitelisted, but were at a Hamming 
        distance > than the treshold chosen. Default: False.
        '''
)


##


# Parse arguments
args = my_parser.parse_args()
path_i = args.input
path_o = args.output
path_spikeins_table = args.path_spikeins_table
cumsum = args.cumsum
n_reads = args.n_reads


########################################################################


# Import code
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d 
import seaborn as sns
import matplotlib.pyplot as plt


# Define helper functions


##


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

    x = df['log10_read_count']
    y = df['log10_n_cells']
    f = interp1d(x, y)
    r_counts_cells = np.corrcoef(x, y)[0,1]
    r_FCs = np.corrcoef(df['FC_from_reads'], df['FC_from_n_cells'])[0,1]

    fig, ax = plt.subplots(figsize=(5.5,5))
    ax.plot(x, y, 'ko')
    sns.regplot(data=df, x='log10_read_count', y='log10_n_cells', scatter=False, ax=ax)
    title = f'Read_counts vs n_cells: r={r_counts_cells:.2f} \n FC from read_counts vs FC from n_cells: r={r_FCs:.2f}'
    ax.set(xlabel='log10 read_count', ylabel='log10 n_cells', title=title)
    ax.spines[['right', 'top']].set_visible(False)
    fig.tight_layout()
    
    return f, fig


##


def filter_GBCs(df_counts, df_spike=None, f=None, path_o=None, use_spike=True, cumsum=75, n_reads=1000):
    """
    Function to filter a read count table with and without spikeins info.
    """
        
    # Calculate log10 read_counts
    df = df_counts.copy()
    
    # W/i
    if use_spike:

        # Remove GBCs falling outside the interpolation range
        min_ = df_spike['log10_read_count'].min()
        max_ = df_spike['log10_read_count'].max()
        df = df.query('log10_read_count >= @min_ and log10_read_count <= @max_')

        # Fit trend
        df['log10_n_cells'] = df['log10_read_count'].map(lambda x: f(x))
        df['n_cells'] = np.round(10 ** df['log10_n_cells'])
        df['cellular_prevalence'] = df['n_cells'] / df['n_cells'].sum()
        df = df.loc[:, ~df.columns.str.contains('log_10|cells')]

    else:
        
        read_count_cumsum = df['read_count'].cumsum()
        t = np.percentile(read_count_cumsum, cumsum)
        gbc_to_retain = read_count_cumsum.loc[lambda x: x<=t].index
        df = df.loc[gbc_to_retain].query('read_count>@n_reads')                             # Filter
        df['cellular_prevalence'] = df['log10_read_count'] / df['log10_read_count'].sum()
        df = df.loc[:, ~df.columns.str.contains('log_10|cells')]

    return df


##


def plot_prevalences(df):
    """
    Plot prevalences.
    """

    # Filter for common
    df_common = df.loc[(df['found_wi'])&(df['found_wo'])]
    x = df_common['cellular_prevalence_wo'].values
    y = df_common['cellular_prevalence_wi'].values
    r = np.corrcoef(x, y)[0,1]

    fig, ax = plt.subplots(figsize=(5.5,5))
    ax.plot(x.cumsum(), 'b.-', label='w/o')
    ax.plot(y.cumsum(), 'r.-', label='w/i')
    ax.legend(loc='upper left', bbox_to_anchor=(1,1), frameon=False, title='Method')
    title = f'Prevalences w/i and w/o spikeins (r={r:.2f})'
    ax.set(xlabel='Ranked GBC', ylabel='Prevance', title=title)
    ax.spines[['right', 'top']].set_visible(False)
    fig.tight_layout()
    
    return fig


##


########################################################################


def main():

    # Read GBC_counts
    df_counts = pd.read_csv(path_i, index_col=0)

    # Filter GBCs
    allowed_status = ['corrected', 'not_whitelisted']
    if args.with_degenerated:
        allowed_status.append('degenerated_not_to_remove')
    df_counts = df_counts.loc[df_counts['status'].isin(allowed_status)]

    # With spikeins
    if path_spikeins_table is not None:

        if os.path.exists(path_spikeins_table):
            print(f'Spikeins detected... Use this info for clonal inference.')
        else:
            raise ValueError(f'Provided path {path_spikeins_table} does not exists!')

    # At least two GBCs are needed to make the interpolation. Better all of them though...
    df_spike = pd.read_csv(path_spikeins_table, index_col=0)
    df_spike.index = df_spike.index.map(lambda x: rev_complement(x))
    assert np.sum(df_spike.index.isin(df_counts.index)) >= 2

    # Fit spikeins trend
    spikeins = df_spike.index[df_spike.index.isin(df_counts.index)]
    df_spike = df_spike.loc[spikeins,:].join(df_counts.loc[spikeins,:])
    df_spike['log10_n_cells'] = np.log10(df_spike['n_cells'])
    df_spike['log10_read_count'] = np.log10(df_spike['read_count'])
    df_spike['FC_from_reads'] = df_spike['read_count'] / df_spike['read_count'].min()
    df_spike['FC_from_n_cells'] = df_spike['n_cells'] / df_spike['n_cells'].min()  # Normalize to min values
    f, fig = fit_trend(df_spike)

    # Save 
    df_spike.to_csv(os.path.join(path_o, 'df_spikeins.csv'))
    fig.savefig(os.path.join(path_o, 'spikeins_fit.png'), dpi=300)

    ##

    # Filter spikeins and plot distributions
    is_spike = df_counts.index.isin(df_spike.index)
    df_counts = df_counts[~is_spike].copy()
    df_counts['log10_read_count'] = np.log10(df_counts['read_count'])

    # W/i
    df_with = filter_GBCs(df_counts, df_spike=df_spike, f=f, path_o=path_o, use_spike=True)

    # W/o
    df_without = filter_GBCs(df_counts, path_o=path_o, use_spike=False, cumsum=cumsum, n_reads=n_reads)

    # Merge info
    df = (
        df_without
        .join(df_with.loc[:,['cellular_prevalence']], rsuffix='_wo', lsuffix='_wi', how='outer')
        .assign(
            found_wi=lambda x: ~x['cellular_prevalence_wi'].isna(),
            found_wo=lambda x: ~x['cellular_prevalence_wo'].isna()
        )
        .sort_values('log10_read_count', ascending=False)
    )
    df.to_csv(os.path.join(path_o, 'clonal_prevalences.csv'))

    # Last vizualization
    fig = plot_prevalences(df)
    fig.savefig(os.path.join(path_o, 'prevalences.png'), dpi=300)

##


########################################################################

# Run
if __name__ == '__main__':
    main()