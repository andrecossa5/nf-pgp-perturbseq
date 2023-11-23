// GET_GBC_ELEMENTs module

nextflow.enable.dsl = 2

//

process GET_GBC_ELEMENTS {

  tag "${sample_name}"
   
  input:
  tuple val(sample_name), path(R1), path(R2), path(filtered)

  output:
  tuple val(sample_name),  path('filtered_R1.fq.gz'), path('filtered_R2.fq.gz'), emit: reads
  tuple val(sample_name),  path('GBC_read_elements.tsv'), emit: elements

  script:
  """
  python ${baseDir}/bin/sc/filter_and_extact_from_GBC_reads.py ${R1} ${R2} ${filtered}
  """

  stub:
  """
  echo ${sample_name} > sample
  touch filtered_R1.fq.gz
  touch filtered_R2.fq.gz
  touch GBC_read_elements.tsv
  """

}

