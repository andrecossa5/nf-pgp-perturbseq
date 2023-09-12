// FASTA_FROM_REF module

nextflow.enable.dsl = 2

//

process FASTA_FROM_REF {

  tag "${sample_name}"
  
  input:
  tuple val(sample_name), val(read_counts)

  output:
  tuple val(sample_name), path("GBC_reference.fa"), emit: fasta

  script:
  """
  awk 'FNR > 1 {print ">"NR-1"\\n"\$1}' ${read_counts} > GBC_reference.fa
  """
  
  stub:
  """
  echo ${sample_name} > sample
  touch GBC_reference.fa
  """

}