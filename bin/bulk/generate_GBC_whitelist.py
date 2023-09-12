import argparse
import pandas as pd
import time
import umi_tools

if __name__ == '__main__':

  parser = argparse.ArgumentParser()

  parser.add_argument(
    '--input',
    type    = str,
    default = '',
    help    = 'Path to input data.'
  )
  parser.add_argument(
    '--method',
    type    = str,
    default = 'directional',
    help    = 'Method to correct barcodes.'
  )
  parser.add_argument(
    '--threshold',
    type    = int,
    default = 1,
    help    = 'Maximum number of nucleotides to correct for.'
  )
  parser.add_argument(
    '--output',
    type    = str,
    default = '',
    help    = 'Path to write whitelist to.'
  )
  args = parser.parse_args()

  print('[' + time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime()) + '] ' + "Loading data from '" + args.input + "'")

  ## read GBCs
  GBC_reads = pd.read_csv(args.input, sep = '\t', names = ['GBC'])

  ## calculate read count for each unique GBC
  GBC_counts = pd.DataFrame(GBC_reads['GBC'].value_counts())
  ## rename column
  GBC_counts.columns = ['read_count']

  ## set method used for clustering GBCs
  clusterer = umi_tools.UMIClusterer(cluster_method = args.method)

  ## convert GBCs to byte type
  GBC_byte = [el.encode('UTF-8') for el in GBC_counts.index.values]

  ## create dict from byte type GBC and read count
  dict = dict(zip(GBC_byte, GBC_counts['read_count']))

  print('[' + time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime()) + '] ' + str(len(dict)) + ' unique GBCs in input data...')
  print('[' + time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime()) + '] ' + "Generating GBC whitelist using the '" + args.method + "' method with a treshold of " + str(args.threshold) + '...')

  ## cluster/group GBCs using the specified data, method, and treshold
  grouped_UMI = clusterer(dict, threshold = args.threshold)

  print('[' + time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime()) + '] ' + str(len(grouped_UMI)) + ' unique GBCs after correction...')
  print('[' + time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime()) + '] ' + "Writing whitelist to '" + args.output + "'")

  ## write list of non-corrected and corrected GBCs to file
  output_file = open(args.output, 'w+')
  for i in range(1, len(grouped_UMI)+1):
    output_file.write(grouped_UMI[i-1][0].decode('utf-8') + '\t' + ','.join([el.decode('utf-8') for el in grouped_UMI[i-1][1:len(grouped_UMI[i-1])]]) + '\n')
  output_file.close()
