// COUNT module

nextflow.enable.dsl = 2

//

// Counts unique (corrected) GBCs reads 
process COUNT {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(GBC_corrected)
  tuple val(sample_name), path(GBC_whitelist) 
  tuple val(sample_name), path(GBC_not_corrected_18bp)

  output:
  tuple val(sample_name), path('read_count_by_GBC_corrected.tsv'), emit: read_counts

  script:
  """
  Rscript \
  ${baseDir}/bin/bulk/count_reads_for_each_GBC.R \
  ${GBC_corrected} \
  ${GBC_whitelist} \
  ${GBC_not_corrected_18bp}
  """

  stub:
  """
  echo ${sample_name} > sample
  touch read_count_by_GBC_corrected.tsv
  """

}

