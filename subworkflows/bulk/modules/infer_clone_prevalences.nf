// INFER_PREVALENCES module

nextflow.enable.dsl = 2

//

// Infere GBC clones frequencies from their read counts using ground truth spikins data
process INFER_PREVALENCES {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(read_count_by_GBC_corrected)

  output:
  tuple val(sample_name), path('clonal_prevalences.csv'), emit: stats_table
  tuple val(sample_name), path('good_GBCs_bulk.txt'), emit: good_GBCs
  tuple val(sample_name), path('spikeins_fit.png'), emit: plot

  script:
  """
  python3 \
  ${baseDir}/bin/bulk/infer_clone_prevalences.py \
  ${params.spikeins_table} \
  ${read_count_by_GBC_corrected}
  """

  stub:
  """
  echo ${sample_name} > sample
  touch clonal_prevalences.csv
  touch good_GBCs_bulk.txt
  touch spikeins_fit.png
  """

}
