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
	default=None,
	help='The path to the output folder of choice. Default: None.'
)
 
 
# Parse arguments
args = my_parser.parse_args()
path_input = args.input
path_output = args.output

# path_input = '/Users/IEO5505/Desktop/prova/prova'
# path_output = '/Users/IEO5505/Desktop/prova'

##


# Import code
import os
import numpy as np
import pandas as pd
import shutil

from Cellula._utils import *


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


##########################################################


# Main
def main():

	# Create folder in path_output
	make_folder(path_output, 'summary', overwrite=True)
	make_folder(os.path.join(path_output, 'summary'), 'prevalences', overwrite=True)
	make_folder(os.path.join(path_output, 'summary'), 'spikeins', overwrite=True)


	folder_d = {}
	for name in os.listdir(path_input):
		path_ = os.path.join(path_input, name)
		if not os.path.islink(path_) and name != '.DS_Store':
			folder_d[name] = path_

	# all_summary
	with open(os.path.join(path_output, 'summary', 'all_summary.txt'), 'wb') as outfile:
		for name in folder_d:
			outfile.write('\n'.encode())
			outfile.write('-------------------------------------'.encode())
			outfile.write('\n'.encode())
			with open(os.path.join(folder_d[name], 'run_summary.txt'), 'rb') as infileobj:
				shutil.copyfileobj(infileobj, outfile)

	# images
	for name in folder_d:
		path_ = os.path.join(path_input, name)
		for x in os.listdir(path_):
			if x == 'spikeins_fit.png':
				shutil.copy(os.path.join(path_, x), os.path.join(path_output, 'summary', 'spikeins', f'{name}.png'))
			elif x == 'prevalences.png':
				shutil.copy(os.path.join(path_, x), os.path.join(path_output, 'summary', 'prevalences', f'{name}.png'))
			elif x == 'df_spikeins.csv':
				shutil.copy(os.path.join(path_, x), os.path.join(path_output, 'summary', 'spikeins', f'{name}.csv'))

	# Prevalences df	
	(	
		concat_df(path_input, folder_d, 'clonal_prevalences')
		.to_csv(os.path.join(path_output, 'summary', 'all_prevalences.csv'))
	)
	# all read_counts df	
	df_correction = concat_df(path_input, folder_d, 'correction_df')
	df_spikeins = concat_df(path_input, folder_d, 'df_spikeins')

	# Write spikeins correction df
	spikeins = df_spikeins.index.unique()
	df_ = (
		df_correction.loc[spikeins]
		.join(df_spikeins.loc[:,['n_cells']])
		.drop_duplicates()
		.sort_values('n_cells', ascending=False)
		.to_csv(os.path.join(path_output, 'summary', 'spikeins_correction.csv'))
	)

##


##########################################################

# Run 
if __name__ == '__main__':
	main()


## Todo:

# # Diagnostic spikeins
# for name in folder_d:
# 	path_ = os.path.join(path_input, name)
# 
# counts = pd.read_csv(os.path.join(path_input, name, 'read_count_by_GBC_corrected.tsv'), index_col=0, sep='\t')
# spikeins = pd.read_csv(os.path.join(path_input, name, 'df_spikeins.csv'), index_col=0).index
# whitelist = (
# 	pd.read_csv(
# 		os.path.join(path_, 'GBC_whitelist.tsv'), 
# 		sep='\t', header=None, names=['GBC', 'similar']
# 	)
# 	.set_index('GBC')
# )
# 
# spikeins_families = whitelist.index[whitelist.index.isin(spikeins)]
# whitelist['similar'] = whitelist['similar'].map(lambda x: x.split(',') if isinstance(x, str) else [])
# 
# 
# df_ = pd.DataFrame({
# 	# 'counts_similar' : whitelist['similar'].map(lambda x: len(x)),
# 	# 'n_spikeins' : whitelist['similar'].map(lambda x: spikeins.isin(pd.Series(x)).sum() ),
# 	'n_corrected' : whitelist['similar'].map(lambda x: pd.Series(x)[pd.Series(x).isin(counts.index)].to_list() ),
# 	'n_uncorrected' : whitelist['similar'].map(lambda x: pd.Series(x)[~pd.Series(x).isin(counts.index)].to_list()  ),
# })
# 
# len(df_.loc['TGCAGTTTTGGTGCTCTA','n_corrected'])
# len(df_.loc['TGCAGTTTTGGTGCTCTA','n_uncorrected'])