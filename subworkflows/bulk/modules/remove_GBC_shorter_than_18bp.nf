// REMOVE_SHORT module

nextflow.enable.dsl = 2

//

// Process
process REMOVE_SHORT {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(GBC_not_corrected)

  output:
  tuple val(sample_name), path('GBC_not_corrected_18bp.tsv.gz'), emit: GBC

  script:
  """
  Rscript \
  ${baseDir}/bin/bulk/remove_GBC_shorter_than_18bp.R \
  ${GBC_not_corrected}
  """

  stub:
  """
  echo ${sample_name} > sample
  touch GBC_not_corrected_18bp.tsv.gz
  """

}