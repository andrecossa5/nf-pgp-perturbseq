// GET_GBC_ELEMENTs module

nextflow.enable.dsl = 2

//

process GET_GBC_ELEMENTS {

  tag "${sample_name}"
   
  input:
  tuple val(sample_name), path(R1), path(R2), path(filtered)

  output:
  tuple val(sample_name),  path('GBC_read_elements.tsv'), emit: elements

  script:
  """
  python ${baseDir}/bin/sc/get_CBC_GBC_UMI.py \
  ${R1} ${R2} ${filtered} ${params.sc_anchor} ${params.sc_anchor_treshold}
  """

  stub:
  """
  echo ${sample_name} > sample
  touch GBC_read_elements.tsv
  """

}

