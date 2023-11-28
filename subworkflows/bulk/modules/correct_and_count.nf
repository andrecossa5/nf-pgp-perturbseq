// CORRECT_AND_COUNT module

nextflow.enable.dsl = 2

//

// Correct GBCs and count their reads
process CORRECT_AND_COUNT {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(GBC)

  output:
  tuple val(sample_name), path('GBC_counts.csv'), emit: counts
  tuple val(sample_name), path('correction_df.csv'), emit: correction_df
  tuple val(sample_name), path('whitelist.csv'), emit: whitelist

  script:
  """
  python \
  ${baseDir}/bin/bulk/correct_and_count.py \
  -i ${GBC} \
  -t ${params.hamming_distance_treshold} \
  --method directional
  """

  stub:
  """
  touch GBC_counts.csv
  touch correction_df.csv
  touch whitelist.csv
  """

}

