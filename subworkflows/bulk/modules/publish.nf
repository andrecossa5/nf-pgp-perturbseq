// publish_bulk module

nextflow.enable.dsl = 2

//

process publish_bulk {

  tag "${sample_name}"

  // Publish
  publishDir "${params.bulk_outdir}/${sample_name}/", mode: 'copy'

  input:
  tuple val(sample_name), path(counts)
  tuple val(sample_name), path(correction_df)
  tuple val(sample_name), path(whitelist)
  tuple val(sample_name), path(stats_table)
  tuple val(sample_name), path(prevalences_plot)
  tuple val(sample_name), path(spikeins_plot)
  tuple val(sample_name), path(df_spikeins)
  tuple val(sample_name), path(run_summary)

  output:
  path counts
  path correction_df
  path whitelist
  path stats_table
  path prevalences_plot
  path spikeins_plot
  path df_spikeins
  path run_summary

  script:
  """
  echo "Moving all necessary files to ${params.bulk_outdir}/${sample_name}/..."
  """

}
