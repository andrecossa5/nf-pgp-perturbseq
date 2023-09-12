// CORRECT module

nextflow.enable.dsl = 2

//

process CORRECT {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(GBC_not_corrected_18bp)
  tuple val(sample_name), path(reformatted_GBC_whitelist)

  output:
  tuple val(sample_name), path('GBC_corrected.tsv.gz'), emit: GBC

  script:
  """
  perl \
  ${baseDir}/bin/bulk/correct_barcodes.pl \
  ${reformatted_GBC_whitelist} \
  ${GBC_not_corrected_18bp} \
  | pigz --fast -p ${task.cpus} \
  > GBC_corrected.tsv.gz
  """

  stub:
  """
  echo ${sample_name} > sample
  touch GBC_corrected.tsv.gz
  """

}