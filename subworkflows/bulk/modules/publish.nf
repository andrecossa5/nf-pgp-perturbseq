// publish_bulk module

nextflow.enable.dsl = 2

//

process publish_bulk {

  tag "${sample_name}"

  // Publish
  publishDir "${params.bulk_outdir}/${sample_name}/", mode: 'copy'

  input:
  tuple val(sample_name), path(GBC_corrected)
  tuple val(sample_name), path(read_count_by_GBC_corrected)
  tuple val(sample_name), path(whitelist)
  tuple val(sample_name), path(clonal_prevalences)
  tuple val(sample_name), path(good_GBCs)
  tuple val(sample_name), path(spikeins_fit)
  tuple val(sample_name), path(run_summary)

  output:
  path GBC_corrected
  path read_count_by_GBC_corrected
  path whitelist
  path clonal_prevalences
  path good_GBCs
  path spikeins_fit
  path run_summary

  script:
  """
  echo "Moving all necessary files to ${params.bulk_outdir}/${sample_name}/..."
  """

}
