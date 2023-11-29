// CORRECT_AND_COUNT module

nextflow.enable.dsl = 2

//

// Correct GBCs and count their reads
process CORRECT_AND_COUNT {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(GBC)

  output:
  tuple val(sample_name), path('GBC_raw_counts.csv.gz'), emit: raw_counts
  tuple val(sample_name), path('GBC_counts_corrected.csv'), emit: corrected_counts
  tuple val(sample_name), path('correction_df.csv'), emit: correction_df

  script:
  """
  python \
  ${baseDir}/bin/bulk/correct_and_count.py \
  -i ${GBC} \
  -t ${params.bulk_graph_clustering_hamming_treshold} \
  --method directional \
  --min_n_reads ${params.bulk_min_n_reads} \
  --min_n_reads ${params.bulk_min_n_reads} \
  --spikeins ${params.spikeins_table}
  """

  stub:
  """
  touch GBC_raw_counts.csv.gz
  touch correction_df.csv
  touch GBC_counts_corrected.csv
  """

}

