// FORMAT_WHITELIST module

nextflow.enable.dsl = 2

//

process FORMAT_WHITELIST {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(GBC_whitelist)

  output:
  tuple val(sample_name), path('GBC_whitelist_proper_format.tsv'), emit: formatted_whitelist

  script:
  """
  Rscript \
  ${baseDir}/bin/bulk/reformat_GBC_whitelist_for_error_correction.R \
  ${GBC_whitelist}
  """

  stub:
  """
  echo ${sample_name} > sample
  touch GBC_whitelist_proper_format.tsv
  """

}