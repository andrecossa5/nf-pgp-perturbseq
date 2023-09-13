// CORRECT_AND_COUNT module

nextflow.enable.dsl = 2

//

// Correct GBCs and count their reads
process CORRECT_AND_COUNT {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(raw_counts)

  output:
  tuple val(sample_name), path('GBC_counts.csv'), emit: counts
  tuple val(sample_name), path('correction_df.csv'), emit: correction_df
  tuple val(sample_name), path('whitelist.csv'), emit: whitelist

  script:
  """
  python3 \
  ${baseDir}/bin/bulk/correct_and_count.py \
  -i ${raw_counts} \
  -t ${params.hamming_distance_treshold} \
  --method directional \
  --min_n_reads ${params.min_n_reads_treshold}
  """

  stub:
  """
  echo ${sample_name} > sample
  touch GBC_counts.csv
  touch correction_df.csv
  touch whitelist.csv
  """

}

