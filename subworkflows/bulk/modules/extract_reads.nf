// EXTRACT_READS module

nextflow.enable.dsl = 2

//

process EXTRACT_READS {

  tag "${sample_name}"

  input:
  tuple val(sample_name), val(in_folder)

  output:
  tuple val(sample_name), path('reads.tsv.gz'), emit: reads

  script:
  """
  zcat ${in_folder}/*.fastq.gz | awk 'NR % 4 == 2' | pigz --fast -p ${task.cpus} > reads.tsv.gz
  """

  stub:
  """
  echo ${sample_name} > sample
  touch reads.tsv.gz
  """
  
}