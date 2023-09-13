#!/usr/bin/python

# Infer clone prevalences script

########################################################################

# Parsing CLI args 
 
# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='GBC_clustering',
    description=
    """
    Script for GBC clustering.
    """
)


# Input
my_parser.add_argument(
    '-i', 
    '--input', 
    type=str,
    default=os.getcwd(),
    help='The path to raw barcodes counts. Default: . .'
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
    help='Maximum Hamming distance for which ... .'
)

# treshold
my_parser.add_argument(
    '--method',
    type=str,
    default='directional',
    help='UMI-tools method. Default: directional.'
)


##


# Parse arguments
args = my_parser.parse_args()
path_i = args.path_input
path_o = args.path_output
method = args.method
threshold = args.threshold


##


# Import code
import numpy as np
import pandas as pd
from umi_tools import UMIClusterer


##


########################################################################

def main():

    # Read and count GBCs
    df = pd.read_csv(path_i, sep='\t', names='GBC')
    counts = df['GBC'].value_counts()
    counts.index = counts.index.map(lambda x : x.encode('UTF-8'))

    # UMI-tools
    clusterer = UMIClusterer(cluster_method=method)
    groups = clusterer(counts.to_dict(), threshold=threshold)
    
    # Extract info and write whitelist 
    with open(path_o, 'w') as f:
        for g in groups:
            head = g[0].decode("utf-8")
            family = ','.join([ x.decode("utf-8") for x in g[1:] ])
            f.write(f'{head}\t{family}\n')


##


########################################################################

# Run
if __name__ == '__main__':
    main()