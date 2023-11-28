#!/usr/bin/python

# Create run summary

########################################################################

# Parsing CLI args 
 
# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='create_run_summary',
    description=
    """
    Create run summary.
    """
)

# indir
my_parser.add_argument(
    '--indir', 
    type=str,
    default=None,
    help='Path to input directory: params.bulk_indir.'
)

# outdir
my_parser.add_argument(
    '--outdir', 
    type=str,
    default=None,
    help='Path to input directory: params.bulk_outdir.'
)

# anchor_sequence
my_parser.add_argument(
    '--anchor_sequence', 
    type=str,
    default=None,
    help='Anchor sequence: params.anchor_sequence.'
)

# read_counts
my_parser.add_argument(
    '--read_counts',
    type=str,
    default=None,
    help='Path to formatted GBC read counts. Default: . .'
)

# Sample
my_parser.add_argument(
    '--sample',
    type=str,
    default=None,
    help='Sample name.'
)

# correction_df
my_parser.add_argument(
    '--correction_df',
    type=str,
    default=None,
    help='Correction_df as generated in correct_and_count.py script.'
)

# stats_table
my_parser.add_argument(
    '--stats_table',
    type=str,
    default=None,
    help='Prevalence df as generated in infer_clone_prevalences.py.'
)

# treshold
my_parser.add_argument(
    '--output',
    type=str,
    default=os.getcwd(),
    help='Ouput folder path. Default: . .'
)


##


# Parse arguments
args = my_parser.parse_args()
indir = args.indir
outdir = args.outdir
anchor_sequence = args.anchor_sequence
path_read_counts = args.read_counts
sample = args.sample
path_correction_df = args.correction_df
stats_table = args.stats_table
path_o = args.output


##


# Import code
import getpass
import datetime
import numpy as np
import pandas as pd


##


########################################################################

def main():

    # Calculate stats

    # Others
    df_ = pd.read_csv(path_read_counts, index_col=0)
    total_GBC_reads = df_['read_count'].sum()
    n_unique_GBCs = df_.shape[0]

    status_counts = df_['status'].value_counts()
    try:
        n_corrected = status_counts['corrected']
        n_not_whitelisted = status_counts['not_whitelisted']
        n_degenerated_removed = status_counts['degenerated_to_remove']
        n_degenerated_not_removed = status_counts['degenerated_not_to_remove']
    except:
        n_corrected = 0
        n_not_whitelisted = 0
        n_degenerated_removed = 0
        n_degenerated_not_removed = 0
        print('Some problem occurred...')
    
    df_ = pd.read_csv(path_correction_df, index_col=0)   
    median_n_degenerated = df_['n_degenerated'].median()
    median_n_reads_added = df_['n_reads_added'].median()

    df_ = pd.read_csv(stats_table, index_col=0) 
    n_unique_wi = np.sum(df_['found_wi'])
    n_unique_wo = np.sum(df_['found_wo'])
    median_pr_wi = df_['cellular_prevalence_wi'].median()
    median_pr_wo = df_['cellular_prevalence_wi'].median()


    ##


    # Write summary
    f = open(os.path.join(path_o, 'run_summary.txt'), 'w')

    # Write stats
    f.write(f'Summary Step 1 (bulk GBC), sample {sample} \n')
    f.write('-------------------------------------')
    f.write('\n')
    f.write('Overview: \n')
    f.write(f'- Date of analysis:               {datetime.datetime.now().strftime("%d-%m-%Y")} \n')
    f.write(f'- User:                           {getpass.getuser()} \n')
    f.write(f'- Working directory:              {"/".join(os.getcwd().split("/")[:-2])} \n')
    f.write('\n')
    f.write(f'Parameters \n')
    f.write(f'--indir:                          {indir} \n')
    f.write(f'--outdir:                         {outdir} \n')
    f.write(f'--anchor_sequence:                {anchor_sequence} \n')
    f.write('\n')
    f.write('Numbers: \n')
    f.write(f'- Reads with GBC (>1 read):       {int(total_GBC_reads)} \n')
    f.write(f'- Unique GBCs:                    {int(n_unique_GBCs)} \n')
    f.write(f'  * corrected:                    {int(n_corrected)} \n')
    f.write(f'  * not_whitelisted:              {int(n_not_whitelisted)} \n')
    f.write(f'  * degenerated, removed:         {int(n_degenerated_removed)} \n')
    f.write(f'  * degenerated, not removed :    {int(n_degenerated_not_removed)} \n')
    f.write(f'- Correction:                          \n')
    f.write(f'  * Median n degenerated GBCs:    {median_n_degenerated:.2f} \n')
    f.write(f'  * Median n added reads:         {median_n_reads_added:.2f} \n')
    f.write(f'- Clonal inference:                    \n')
    f.write(f'  *  w/i spikeins:                {int(n_unique_wi)} unique clones, median prevalence {median_pr_wi:.2f}  \n')
    f.write(f'  *  w/o spikeins:                {int(n_unique_wo)} unique clones, median prevalence {median_pr_wo:.2f}  \n')
    f.write('\n')

    f.close()
 

##


########################################################################

# Run
if __name__ == '__main__':
    main()