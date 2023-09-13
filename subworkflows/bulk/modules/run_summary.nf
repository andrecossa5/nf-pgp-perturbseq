// generate_run_summary_bulk module

nextflow.enable.dsl = 2

//

// Process
process generate_run_summary_bulk {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(reads)
  tuple val(sample_name), path(read_counts)
  tuple val(sample_name), path(correction_df)
  tuple val(sample_name), path(stats_table)

  output:
  tuple val(sample_name), path('run_summary.txt'), emit: summary

  script:
  """
  python3 \
  ${baseDir}/bin/bulk/create_run_summary.py \
  --reads ${reads} \
  --indir ${params.indir} \
  --outdir ${params.outdir} \
  --anchor_sequence ${params.anchor_sequence} \
  --sample ${sample} \
  --read_counts ${read_counts} \
  --correction_df ${correction_df} \
  --stats_table ${stats_table} 
  """

  stub:
  """
  echo ${sample_name} > sample
  touch run_summary.txt
  """

}