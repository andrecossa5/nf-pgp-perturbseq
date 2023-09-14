#!/usr/bin/python

# Correct and count genomic GBC

########################################################################

# Parsing CLI args 
 
# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='correct_and_count',
    description=
    """
    Script for GBC filtering, correction and counting.
    """
)


# Input
my_parser.add_argument(
    '-i', 
    '--input', 
    type=str,
    default=None,
    help='''
        The path to extracted GBCs file. Must be a .csv file, 
        with one column. Default: . .
        '''
)

# Output
my_parser.add_argument(
    '-o', 
    '--output', 
    type=str,
    default=os.getcwd(),
    help='Output folder path. Default: . .'
)

# treshold
my_parser.add_argument(
    '-t', 
    '--threshold',
    type=int,
    default=1,
    help='Maximum Hamming distance for which sequences can be clustered together.'
)

# treshold
my_parser.add_argument(
    '--method',
    type=str,
    default='directional',
    help='UMI-tools method. Default: directional.'
)

# treshold
my_parser.add_argument(
    '--min_n_reads',
    type=int,
    default=100,
    help='min_n_reads that a GBCs must have to be considered in the final counts. Default: 100.'
)


##


# Parse arguments
args = my_parser.parse_args()
path_i = args.input
path_o = args.output
method = args.method
threshold = args.threshold
min_n_reads = args.min_n_reads


# path_i = '/Users/IEO5505/Desktop/prova/GBC_not_corrected.tsv.gz'
# path_o = '/Users/IEO5505/Desktop/prova'
# method = 'directional'
# threshold = 3
# min_n_reads = 100


##


# Import code
import numpy as np
import pandas as pd
from umi_tools import UMIClusterer
from scipy.spatial.distance import hamming


##


# Helpers
def process_one_GBC(df, threshold):

    n_degenerated = df.shape[0]
    test = df['hamming']<=threshold
    fract_below_treshold = test.sum() / n_degenerated
    n_before = df['n_reads_correct'].values[0]
    n_reads_after_correction = n_before + df.loc[test, 'n_reads_degenerated'].sum()
    diff = n_reads_after_correction - n_before

    d = {
        'n_degenerated' : n_degenerated,
        'fract_below_treshold' : fract_below_treshold,
        'n_reads_before_correction' : n_before,
        'n_reads_after_correction' : n_reads_after_correction,
        'n_reads_added' : diff
    }

    return pd.Series(d)


##


########################################################################

def main():

    # Remove short GBCs (len == 18bp), count and save all raw_counts
    GBCs = pd.read_csv(path_i, header=None)[0]
    is_18bp = GBCs.map(lambda x: len(x) == 18)
    GBCs = GBCs.loc[is_18bp].copy()
    GBC_counts = GBCs.value_counts().astype(np.int32)
    del GBCs    # Save memory

    # Generate GBC clusters, and reformat
    GBC_counts.index = GBC_counts.index.map(lambda x : x.encode('UTF-8'))
    clusterer = UMIClusterer(cluster_method=method)
    groups = clusterer(GBC_counts.to_dict(), threshold=threshold)

    # To long whitelist df
    df_l = [
        pd.DataFrame(
            [ ( g[0].decode("utf-8"), x.decode("utf-8") ) for x in g[1:] ],
            columns=['correct', 'degenerated']
        )
        for g in groups
    ]
    del groups                  # Save memory
    df = pd.concat(df_l)
    del df_l                    # Save memory

    # Calculate hammings
    df['hamming'] = [
        hamming(np.array(list(x)), np.array(list(y))) * 18 \
        for x, y in zip(df['correct'], df['degenerated'])
    ]

    # Add counts info to whitelist_df
    GBC_counts.index = GBC_counts.index.map(lambda x : x.decode('UTF-8')) # Restore encoding
    df['n_reads_correct'] = GBC_counts.loc[df['correct']].values
    df['n_reads_degenerated'] = GBC_counts.loc[df['degenerated']].values

    # Filter out correct-degenerated pairs with too low or too similar read_counts:
    # This will not be considered for correction.
    df = df.query('n_reads_correct>=@min_n_reads and n_reads_correct>=@min_n_reads*n_reads_degenerated')

    # To the remaining correct GBCs, add counts from the corresponding degenerated barcodes 
    # NB: not from all of them, but only the ones at <= <treshold> hamming distance
    correction_df = df.groupby('correct').apply(lambda x: process_one_GBC(x, threshold))
    GBC_counts.loc[correction_df.index] = correction_df['n_reads_after_correction']

    # Annotate GBC_counts
    GBC_counts = GBC_counts.to_frame('read_count')
    L = [
        GBC_counts.index.isin(correction_df.index),
        GBC_counts.index.isin(df['correct']) & ~GBC_counts.index.isin(correction_df.index),
        GBC_counts.index.isin(df.query('hamming<=@threshold')['degenerated']),
        GBC_counts.index.isin(df.query('hamming>@threshold')['degenerated'])
    ]
    values = [ 'corrected', 'discarded', 'degenerated_to_remove', 'degenerated_not_to_remove' ]
    GBC_counts['status'] = np.select(L, values, default='not_whitelisted')

    # Save
    GBC_counts.to_csv(os.path.join(path_o, 'GBC_counts.csv'))
    correction_df.to_csv(os.path.join(path_o, 'correction_df.csv'))
    df.to_csv(os.path.join(path_o, 'whitelist.csv'))


##


########################################################################

# Run
if __name__ == '__main__':
    main()