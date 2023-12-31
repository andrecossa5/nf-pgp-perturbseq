#!/usr/bin/python

# Create summary output script 

########################################################################

# Parsing CLI args 

# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
	prog='collapse_outputs',
	description=
	"""
	Script for aggregating outputs nf-perturbseq pipeline, bulk.
	"""
)

# Add arguments
 
# Input
my_parser.add_argument(
	'-i', 
	'--input', 
	type=str,
	default=None,
	help='The path to results/perturb_bulk folder. Default: None.'
)
 
# Ouput
my_parser.add_argument(
	'-o', 
	'--output', 
	type=str,
	default=os.getcwd(),
	help='''
		The path to the output folder of choice. Default: os.getcwd().
		'''
)
 
 
# Parse arguments
args = my_parser.parse_args()
path_input = args.input
path_output = args.output

##

# Import code
import numpy as np
import pandas as pd
import shutil
from scipy.stats import spearmanr, pearsonr
from plotting_utils._utils import *


##


# Helpers
def concat_df(path_input, folder_d, string_pattern):

	L = []
	for name in folder_d:
		path_ = os.path.join(path_input, name)
		df_ = pd.read_csv(os.path.join(path_, f'{string_pattern}.csv'), index_col=0)
		L.append(df_.assign(sample=name))
	df = pd.concat(L)

	return df 


##


def get_folder_d(path_input):
		
	folder_d = {}
	for name in os.listdir(path_input):
		if name != '.DS_Store' and name != 'summary':
			folder_d[name] = os.path.join(path_input, name)

	return folder_d


##


def write_all_summary(path_output, folder_d):
	"""
	Summary of all runs. General infos.
	"""
	with open(os.path.join(path_output, 'summary', 'all_summary.txt'), 'wb') as outfile:
		for name in folder_d:
			outfile.write('\n'.encode())
			outfile.write('-------------------------------------'.encode())
			outfile.write('\n'.encode())
			with open(os.path.join(folder_d[name], 'run_summary.txt'), 'rb') as infileobj:
				shutil.copyfileobj(infileobj, outfile)


##


def write_correction_summary(path_input, path_output, folder_d):
	"""
	Write summary of every run correction procedure.
	"""
	d = {}
	for name in folder_d:
		d_ = {}
		df_all = pd.read_csv(os.path.join(path_input, name, 'GBC_counts.csv'), index_col=0)
		df_good = pd.read_csv(os.path.join(path_input, name, 'clonal_prevalences.csv'), index_col=0)
		df_corr = pd.read_csv(os.path.join(path_input, name, 'correction_df.csv'), index_col=0)
		d_['n'] = df_all.shape[0]
		d_['JI'] = np.sum(df_good['found_wi'] & df_good['found_wo']) / df_good.shape[0]
		d_['perc_corrected'] = df_all.query('status == "corrected"').shape[0] / df_all.shape[0]
		d_['perc_added_reads'] = np.median(df_corr['n_reads_added'] / df_corr['n_reads_before_correction'])
		df_common = df_good.loc[lambda x: x['found_wi'] & x['found_wo']]
		if df_common.shape[0]>10:
			d_['pearson'] = pearsonr(df_common['cellular_prevalence_wi'], df_common['cellular_prevalence_wo'])[0]
			d_['spearman'] = spearmanr(df_common['cellular_prevalence_wi'], df_common['cellular_prevalence_wo'])[0]
		d[name] = d_
	df = pd.DataFrame(d).T

	# Write trace
	with open(os.path.join(path_output, 'summary', 'correction_summary.txt'), 'w') as f:
		f.write('---------------------------------------------- \n')
		f.write(f'# Summary correction procedure across samples \n')
		f.write('---------------------------------------------- \n')
		f.write('\n')
		f.write('- All numbers presented across sample (median +- std): \n')
		f.write('\n')
		f.write(f' * n unique barcodes (all barcodes): {df["n"].median():.2f} (+-{df["n"].std():.2f}) \n')
		f.write(f' * J.I. clones inferred w/i and w/o spikeins: {df["JI"].median():.2f} (+-{df["JI"].std():.2f}) \n')
		f.write(f' * % of corrected barcodes: {df["perc_corrected"].median():.2f} (+-{df["perc_corrected"].std():.2f}) \n')
		f.write(f' * % reads added at correction: {df["perc_added_reads"].median():.2f} (+-{df["perc_added_reads"].std():.2f}) \n')
		f.write('\n')
		f.write('- Cellular prevalences correlation (i.e., common clones detected w/i and w/o spikeins): \n')
		try:
			f.write(f' * Pearson\'s correlation: {df["pearson"].median():.2f} (+-{df["pearson"].std():.2f}) \n')
			f.write(f' * Spearman\'s correlation: {df["spearman"].median():.2f} (+-{df["spearman"].std():.2f}) \n')
		except:
			f.write(f' * Pearson\'s correlations cannot be calculated...')
		f.write('\n')


##


##########################################################


# Main
def main():

	# Create folders in path_output
	make_folder(path_output, 'summary', overwrite=True)
	make_folder(os.path.join(path_output, 'summary'), 'prevalences', overwrite=True)
	make_folder(os.path.join(path_output, 'summary'), 'spikeins', overwrite=True)

	folder_d = get_folder_d(path_input)

	# Copy necessary files from each run
	for name in folder_d:
		path_ = os.path.join(path_input, name)
		for x in os.listdir(path_):
			if x == 'spikeins_fit.png':
				shutil.copy(os.path.join(path_, x), 
					os.path.join(path_output, 'summary', 'spikeins', f'{name}.png'))
			elif x == 'final_prevalences.png':
				shutil.copy(os.path.join(path_, x), 
					os.path.join(path_output, 'summary', 'prevalences', f'{name}.png'))
			elif x == 'df_spikeins.csv':
				shutil.copy(os.path.join(path_, x),
					os.path.join(path_output, 'summary', 'spikeins', f'{name}.csv'))
			elif x == 'GBC_1read_distributions.png':
				shutil.copy(os.path.join(path_, x),
					os.path.join(path_output, 'summary', 'prevalences', f'GBC_1read_distributions_{name}.png'))
			elif x == 'GBC_morereads_distributions.png':
				shutil.copy(os.path.join(path_, x),
					os.path.join(path_output, 'summary', 'prevalences', f'GBC_morereads_distributions_{name}.png'))

	# Prevalences df	
	df_prevalences = concat_df(path_input, folder_d, 'clonal_prevalences')
	df_prevalences.to_csv(os.path.join(path_output, 'summary', 'all_prevalences.csv'))

	# all read_counts df	
	df_correction = concat_df(path_input, folder_d, 'correction_df')
	df_correction.to_csv(os.path.join(path_output, 'summary', 'df_correction.csv'))
	df_spikeins = concat_df(path_input, folder_d, 'df_spikeins')

	# Write spikeins correction df
	df_ = (
		df_correction.loc[df_correction.index.isin(df_spikeins.index)]
		.join(df_spikeins.loc[:,['n_cells']])
		.drop_duplicates()
		.sort_values('n_cells', ascending=False)
	)
	df_.to_csv(os.path.join(path_output, 'summary', 'spikeins_correction.csv'))

	# Summaries
	write_all_summary(path_output, folder_d)
	write_correction_summary(path_input, path_output, folder_d)


	##


##########################################################

# Run 
if __name__ == '__main__':
	main()
