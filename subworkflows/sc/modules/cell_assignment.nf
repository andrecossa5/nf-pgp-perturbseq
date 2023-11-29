// CELL_ASSIGNMENT module

nextflow.enable.dsl = 2

//

process CELL_ASSIGNMENT {

  tag "${sample_name}"
  
  input:
  tuple val(sample_name), path(elements)

  output:
  tuple val(sample_name), path("CBC_GBC_combos.tsv.gz"), emit: CBC_GBC_combos
  tuple val(sample_name), path("clones_summary_table.csv"), emit: clones_summary
  tuple val(sample_name), path("cells_summary_table.csv"), emit: cells_summary
  tuple val(sample_name), path("CBC_GBC_combo_status.png"), emit: plot

  script: 
  """
  python ${baseDir}/bin/sc/cell_assignment.py \
  --sample ${sample_name} \
  --path_bulk ${params.bulk_outdir}/summary/bulk_GBC_reference.csv \ 
  --path_sc ${elements} \
  --method ${params.cell_assignment_method} \
  --ncores ${task.cpus} \
  --bulk_sc_treshold ${params.bulk_sc_treshold} \
  --umi_treshold ${params.umi_treshold} \
  --read_treshold ${params.read_treshold} \
  --coverage_treshold ${params.coverage_treshold}
  """

  stub:
  """
  echo ${sample_name} > sample
  touch CBC_GBC_combos.tsv.gz
  touch clones_summary_table.csv
  touch cells_summary_table.csv
  touch CBC_GBC_combo_status.png
  """

}