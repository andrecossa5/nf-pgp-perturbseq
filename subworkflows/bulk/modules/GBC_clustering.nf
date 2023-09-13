// WHITELIST module

nextflow.enable.dsl = 2

//

// Process
process WHITELIST {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(GBC_not_corrected_18bp)

  output:
  tuple val(sample_name), path("GBC_whitelist.tsv"), emit: whitelist

  script:
  """
  python3 \
  ${baseDir}/bin/bulk/GBC_clustering.py \
  --i ${GBC_not_corrected_18bp} \
  --method cluster \
  --threshold 1
  """

  stub:
  """
  echo ${sample_name} > sample
  touch GBC_whitelist.tsv
  """

}

