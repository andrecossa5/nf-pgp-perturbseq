// FASTA_FROM_REF module

nextflow.enable.dsl = 2

//

process FASTA_FROM_REF {

  tag "${sample_name}"
   
  input:
  tuple val(sample_name), val(in_folder)

  output:
  tuple val(sample_name), path("GBC_reference.fa"), emit: fasta

  script:
  """
  python3 fasta_from_ref.py -i ${params.bulk_outdir}/${sample_name}
  """
  
  stub:
  """
  echo ${sample_name} > sample
  touch GBC_reference.fa
  """

}