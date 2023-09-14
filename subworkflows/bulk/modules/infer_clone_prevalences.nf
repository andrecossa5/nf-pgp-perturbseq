// INFER_PREVALENCES module

nextflow.enable.dsl = 2

//

// Infere GBC clones frequencies from their read counts using ground truth spikins data
process INFER_PREVALENCES {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(counts)

  output:
  tuple val(sample_name), path('df_spikeins.csv'), emit: df_spikeins
  tuple val(sample_name), path('spikeins_fit.png'), emit: spikeins_plot
  tuple val(sample_name), path('clonal_prevalences.csv'), emit: stats_table
  tuple val(sample_name), path('prevalences.png'), emit: prevalences_plot

  script:
  """
  python3 \
  ${baseDir}/bin/bulk/infer_clone_prevalences.py \
  -i ${counts}
  --path_spikeins_table ${params.spikeins_table}
  """

  stub:
  """
  echo ${sample_name} > sample
  touch clonal_prevalences.csv
  touch df_spikeins.csv
  touch spikeins_fit.png
  touch prevalences.png
  """

}
