// FIND_GBC module

nextflow.enable.dsl = 2
params.bin_dir = "../bin/"

//

// Process
process FIND_GBC {

  tag "${sample_name}"

  input:
  path(search_patterns)
  tuple val(sample_name), path(reads)

  output:
  tuple val(sample_name), path('GBC_not_corrected.tsv.gz'), emit: GBC

  script:
  """
  zcat ${reads} \
  | egrep -f ${search_patterns} -o \
  | awk '{print substr(\$0, 23, 18);}' \
  | pigz --fast -p ${task.cpus} \
  > GBC_not_corrected.tsv.gz
  """

  stub:
  """
  echo ${sample_name} > sample
  touch GBC_not_corrected.tsv.gz
  """

}